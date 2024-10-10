#pragma once

#include <chrono>
#include <thread>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <valik/shared.hpp>
#include <valik/search/local_prefilter.hpp>
#include <valik/search/query_record.hpp>
#include <valik/search/sync_out.hpp>
#include <utilities/cart_queue.hpp>
#include <utilities/threshold/basics.hpp>

#include <raptor/threshold/threshold.hpp>

namespace valik::app
{

/**
 * @brief Create parallel prefiltering jobs.
*/
template <typename query_t, typename index_t>
inline void prefilter_queries_parallel(index_t const & index,
                                       search_arguments const & arguments,
                                       std::vector<query_t> const & records,
                                       raptor::threshold::threshold const & thresholder,
                                       cart_queue<query_t> & queue)
{
    if (records.empty())
        return;

    std::vector<std::jthread> tasks;
    size_t const num_records = records.size();
    size_t const records_per_thread = num_records / arguments.threads;

    sync_out verbose_out(arguments.disabledQueriesFile);
    for (size_t i = 0; i < arguments.threads; ++i)
    {
        size_t const start = records_per_thread * i;
        size_t const end = i == (unsigned) (arguments.threads - 1) ? num_records : records_per_thread * (i + 1);

        std::span<query_t const> records_slice{&records[start], &records[end]};

        auto prefilter_cb = [&queue,&arguments,&verbose_out,&index](query_t const & record, 
                                                                  std::unordered_set<size_t> const & bin_hits)
        {
            auto & ibf = index.ibf();
            if (bin_hits.size() > std::max((size_t) 4, (size_t) std::round(ibf.bin_count() / 2.0)))
            {
                if (arguments.verbose)
                    verbose_out.write_warning(record, bin_hits.size());
                if (arguments.keep_best_repeats)    // keep hits for bins with the highest entropy
                {
                    auto const & entropy_ranking = index.entropy_ranking();
                    size_t inserted_bins{0};
                    for (auto bin : entropy_ranking)
                    {
                        if (bin_hits.count(bin) > 0)
                        {
                            queue.insert(bin, record);
                            inserted_bins++;
                        }
                        if (inserted_bins < std::max((size_t) 4, (size_t) std::round(ibf.bin_count() / 10.0)))
                            return;
                    }
                }
                else if (arguments.keep_all_repeats)
                {
                    for (auto & bin : bin_hits)
                    {
                        queue.insert(bin, record);
                    }
                }
                return;
            }
            
            for (auto const bin : bin_hits)
            {
                queue.insert(bin, record);
            }
        };

        // The following calls `local_prefilter(records, ibf, arguments, threshold)` on a thread.
        tasks.emplace_back([=, &index, &arguments, &thresholder]()
        {
            local_prefilter(records_slice, index.ibf(), arguments, thresholder, prefilter_cb);
        });
    }
}

/**
 * @brief Fill the shopping cart queue with all matches. Search each query record in each reference bin without prefiltering.
*/
template <typename query_t>
inline void search_all_parallel(size_t const ref_seg_count,
                                search_arguments const & arguments,
                                std::vector<query_t> const & records,
                                cart_queue<query_t> & queue)
{
    if (records.empty())
        return;

    std::vector<std::jthread> tasks;
    size_t const num_records = records.size();
    size_t const records_per_thread = num_records / arguments.threads;

    for (size_t i = 0; i < arguments.threads; ++i)
    {
        size_t const start = records_per_thread * i;
        size_t const end = i == (unsigned) (arguments.threads - 1) ? num_records : records_per_thread * (i + 1);

        std::span<query_t const> records_slice{&records[start], &records[end]};

        auto all_cb = [=,&queue,&arguments](query_t const& record)
        {
            for (size_t bin{0}; bin < ref_seg_count; bin++)
            {
                queue.insert(bin, record);
            }
        };

        tasks.emplace_back([=]()
        {
            for (query_t const & record : records_slice)
                all_cb(record);
        });
    }
}

/**
 * @brief Create a queue of bin indices to send for distributed search. Search each query record in each reference bin without prefiltering.
 *        Query I/O from disk in external process. 
*/
template <typename query_t>
inline void fill_queue_with_bin_ids(size_t const ref_seg_count,
                                    search_arguments const & arguments,
                                    cart_queue<query_t> & queue)
{
    std::vector<std::jthread> tasks;
    size_t const bins_per_thread = ref_seg_count / arguments.threads;

    query_t dummy_record{};
    for (size_t i = 0; i < arguments.threads; ++i)
    {
        size_t const start = bins_per_thread * i;
        size_t const end = i == (unsigned) (arguments.threads - 1) ? ref_seg_count : bins_per_thread * (i + 1);

        auto bin_cb = [=,&queue](size_t const bin)
        {
            queue.insert(bin, dummy_record);
        };

        tasks.emplace_back([=]()
        {
            for (size_t bin = start; bin < end; bin++)
                bin_cb(bin);
        });
    }
}

} // namespace valik::app
