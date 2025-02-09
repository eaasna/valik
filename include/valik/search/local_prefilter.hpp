#pragma once

#include <span>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/threshold/threshold.hpp>

#include <valik/search/query_record.hpp>
#include <valik/shared.hpp> // search_arguments

namespace valik
{

/**
 * @brief Function that samples patterns on a query.
 *
 * @param read_len Length of query.
 * @param pattern_size Length of pattern.
 * @param query_every Every nth potential match is considered.
 * @param callback Functor that corrects the threshold based on matching k-mer counts.
 * @return Lower quartile of threshold correction.
 */
template <typename functor_t>
constexpr double sample_begin_positions(size_t const read_len, uint64_t const pattern_size, uint8_t const query_every, functor_t && callback)
{
    assert(read_len >= pattern_size);

    size_t first_pos{pattern_size};    
    if (read_len < query_every + pattern_size)
        first_pos = 0;    // start from beginning when short query

    size_t corrected_pattern_count{0u};
    double total_correction{0};
    for (size_t pos = first_pos; pos <= read_len - pattern_size; pos = pos + query_every * pattern_size * 2)
    {
        auto correction = callback(pos);
        if (correction > 0)
        {
            corrected_pattern_count++;
            total_correction += correction;
        }
    }

    return 1 + total_correction / (double) std::max(corrected_pattern_count, (size_t) 1);
}


/**
 * @brief Function that finds the begin positions of all pattern of a query.
 *
 * E.g if read_len = 150
 *        pattern_size = 50
 *        query_every = 30
 *
 *        begin_positions = {0, 30, 60, 90, 100}
 *
 * @param read_len Length of query.
 * @param pattern_size Length of pattern.
 * @param query_every Every nth potential match is considered.
 * @param callback Functor that calls make_pattern_bounds and find_pattern_bins on each begin position.
 */
template <typename functor_t>
constexpr void pattern_begin_positions(size_t const read_len, uint64_t const pattern_size, uint8_t const query_every, functor_t && callback)
{
    assert(read_len >= pattern_size);

    size_t last_begin{0u};
    for (size_t pos = 0; pos <= read_len - pattern_size; pos = pos + query_every)
    {
        callback(pos);
        last_begin = pos;
    }

    if (last_begin < read_len - pattern_size)
    {
        // last pattern might have a smaller overlap to make sure the end of the record is covered
        callback(read_len - pattern_size);
    }
}

/**
 * @brief Struct representing a half open interval [begin, end) of the minimiser vector corresponding to a particular pattern.
 */
struct pattern_bounds
{
    size_t begin_position;
    size_t end_position;
    size_t threshold;

    size_t minimiser_count() const
    {
        return end_position - begin_position;
    }
};

/**
 * @brief Function that matches a pattern to the corresponding range of minimisers in the minimisers vector and finds the threshold based on
 *        the number of minimisers in that range. For winnowing minimisers patterns of the same length can contain a varying number of minimisers.
 *
 * @param begin Begin position of the pattern on the query.
 * @param arguments Command line arguments.
 * @param window_span_begin Vector containing for each minimiser the beginning of the first window to which it is the minimiser.
 * @param thresholder Threshold for the number of shared k-mers to constitute a likely local match.
 * @return pattern_bounds The interval of minimisers corresponding to the pattern.
 */
template <typename span_vec_t>
pattern_bounds make_pattern_bounds(size_t const & begin,
                                   search_arguments const & arguments,
                                   span_vec_t const & window_span_begin,
                                   raptor::threshold::threshold const & thresholder)
{
    assert(window_span_begin.size() >= 1);
    assert(window_span_begin[0] == 0);

    auto pattern = pattern_bounds{};

    // upper bound returns the first element of window_span_begin that is > begin
    // the bound element is the second minimiser of the pattern
    auto begin_it = std::upper_bound(window_span_begin.begin(), window_span_begin.end(), begin);
    assert(begin_it != window_span_begin.end());
    pattern.begin_position = begin_it - window_span_begin.begin() - 1;  // -1 to find the first minimiser of the pattern

    size_t last_window_of_pattern = begin + arguments.pattern_size - arguments.window_size;
    // the bound element is the first minimiser after the pattern
    auto end_it = std::upper_bound(window_span_begin.begin(), window_span_begin.end(), last_window_of_pattern);
    assert(end_it != window_span_begin.begin());
    pattern.end_position = end_it - window_span_begin.begin();

    pattern.threshold = thresholder.get(pattern.minimiser_count());
    return pattern;
}

/**
 * @brief Function that for a single pattern counts matching k-mers and corrects the threshold to avoid too many spuriously matching bins.
 *
 * @param pattern Slice of a query record that is being considered.
 * @param bin_count Number of bins in the IBF.
 * @param counting_table Rows: minimisers of the query. Columns: bins of the IBF.
 * @return Threshold correction that avoids too many spurious matches.
 */
template <typename binning_bitvector_t>
double find_dynamic_threshold_correction(pattern_bounds const & pattern,
                                          size_t const & bin_count,
                                          binning_bitvector_t const & counting_table)
{
    // counting vector for the current pattern
    seqan3::counting_vector<uint8_t> total_counts(bin_count, 0);

    for (size_t i = pattern.begin_position; i < pattern.end_position; i++)
        total_counts += counting_table[i];

    std::unordered_set<size_t> pattern_hits;

    bool max_threshold{false};
    uint8_t correction_count{0};
    while (true)
    {
        for (size_t current_bin = 0; current_bin < total_counts.size(); current_bin++)
        {
            auto &&count = total_counts[current_bin];
            if (count >= (pattern.threshold + correction_count))
            {
                pattern_hits.insert(current_bin);
            }
        }
        if ((pattern.threshold + correction_count) >= pattern.minimiser_count())
            max_threshold = true;
        if (pattern_hits.size() < std::max((size_t) 4, (size_t) std::round(bin_count / 4.0)) || 
            max_threshold)
            break;
        else
        {
            pattern_hits.clear();
            // increase threshold in 10% increments or by at least 1 to find lowest threshold that is not ubiquitous
            correction_count += std::max((size_t) 1, (size_t) std::round(pattern.threshold * 0.1 * correction_count));
        }
    }
    
    return (double) correction_count / (double) pattern.threshold;
}


/**
 * @brief Function that for a single pattern counts matching k-mers and returns bins that exceed the threshold.
 *
 * @param pattern Slice of a query record that is being considered.
 * @param correction Threshold correction determined from a sample of patterns.
 * @param bin_count Number of bins in the IBF.
 * @param counting_table Rows: minimisers of the query. Columns: bins of the IBF.
 * @param sequence_hits Bins that likely contain a match for the pattern (IN-OUT parameter).
 */
template <typename binning_bitvector_t>
void find_pattern_bins(pattern_bounds const & pattern,
                       double const & correction_coef,
                       size_t const & bin_count,
                       binning_bitvector_t const & counting_table,
                       std::unordered_set<size_t> & sequence_hits)
{
    // counting vector for the current pattern
    seqan3::counting_vector<uint8_t> total_counts(bin_count, 0);

    for (size_t i = pattern.begin_position; i < pattern.end_position; i++)
        total_counts += counting_table[i];

    for (size_t current_bin = 0; current_bin < total_counts.size(); current_bin++)
    {
        auto &&count = total_counts[current_bin];
        if (count >= (pattern.threshold * correction_coef))
        {
            // the result is a union of results from all patterns of a read
            sequence_hits.insert(current_bin);
        }
    }
}

/**
 * @brief Function that queries the IBF for local matches in a batch of records.
 *
 * @param records Query records.
 * @param ibf Interleaved Bloom Filter of the reference database.
 * @param arguments Command line arguments.
 *                  arguments.pattern_size and arguments.error_rate define the minimum length and maximum error rate of a local match respectively.
 *                  arguments.query_every defines how many match locations are considered per record.
 * @param thresholder Threshold for the number of shared k-mers to constitute a likely local match.
 * @param result_cb Lambda that inserts the prefiltering results (record-bin pairs) into the shopping carts.
 */
template <seqan3::data_layout ibf_data_layout, typename result_cb_t, typename query_t>
void local_prefilter(
    std::span<query_t const> const & records,
    seqan3::interleaved_bloom_filter<ibf_data_layout> const & ibf,
    search_arguments const & arguments,
    raptor::threshold::threshold const & thresholder,
    result_cb_t result_cb)
{
    // concurrent invocations of the membership agent are not thread safe
    // agent has to be created for each thread
    auto agent = ibf.membership_agent();
    size_t const bin_count = ibf.bin_count();

    // vector holding all the minimisers and their starting position for the read
    std::vector<std::tuple<uint64_t, size_t>> minimiser;

    auto minimiser_hash_adaptor = seqan3::views::minimiser_hash(
        arguments.shape,
        seqan3::window_size{arguments.window_size},
        seqan3::seed{adjust_seed(arguments.shape_weight)});

    for (query_t const & record : records)
    {
        std::vector<seqan3::dna4> const & seq = record.sequence;

        // sequence can't contain local match if it's shorter than pattern length
        if (seq.size() < arguments.pattern_size)
            continue;

        // basically: minimiser = seq | minimiser_hash_adaptor | seqan3::views::to<decltype(minimiser)>;
        {
            auto const minimiser_hash = minimiser_hash_adaptor(seq);
            auto it = minimiser_hash.begin();
            auto const sentinel = minimiser_hash.end();
            auto const hash_begin = it.base();
            for (; it != sentinel; ++it)
                minimiser.emplace_back(*it, it.base() - hash_begin);
        }

        //-----------------------------
        //
        // Table of counting vectors newly created for each read
        // rows: each minimiser of read
        // columns: each bin of IBF
        //
        //-----------------------------
        using binning_bitvector_t = typename std::remove_cvref_t<decltype(ibf)>::membership_agent_type::binning_bitvector;
        std::vector<binning_bitvector_t> counting_table(minimiser.size(),
                                                        binning_bitvector_t(bin_count));

        // the beginning of the first window this minimiser is in
        std::vector<size_t> window_span_begin(minimiser.size(), 0);

        for (size_t i{0}; i < minimiser.size(); i++)
        {
            const auto &[min, start_pos] = minimiser[i];
            window_span_begin[i] = start_pos;
            counting_table[i].raw_data() |= agent.bulk_contains(min).raw_data();
        }

        minimiser.clear();

        std::unordered_set<size_t> sequence_hits{};
        double threshold_correction{1};
        if (!arguments.static_threshold)
        {
            threshold_correction = sample_begin_positions(seq.size(), arguments.pattern_size, arguments.query_every, [&](size_t const begin) -> double
            {
                pattern_bounds const pattern = make_pattern_bounds(begin, arguments, window_span_begin, thresholder);
                return find_dynamic_threshold_correction(pattern, bin_count, counting_table);
            });
        }

        pattern_begin_positions(seq.size(), arguments.pattern_size, arguments.query_every, [&](size_t const begin)
        {
            pattern_bounds const pattern = make_pattern_bounds(begin, arguments, window_span_begin, thresholder);
            find_pattern_bins(pattern, threshold_correction, bin_count, counting_table, sequence_hits);
        });

        result_cb(record, sequence_hits);
    }
}

} // namespace valik
