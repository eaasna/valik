#pragma once

#include <vector>
#include <map>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/std/span>

#include <raptor/threshold/threshold.hpp>

#include <valik/search/query_record.hpp>
#include <valik/search/query_result.hpp>
#include <valik/shared.hpp> // search_arguments

namespace valik
{

//-----------------------------
//
// Reports all pattern begin positions in read.
//
// If read_len = 150
//   pattern_size = 50
//   overlap = 20
//
//   begin_vector = {0, 30, 60, 90, 100}
//
//-----------------------------
template <typename functor_t>
constexpr void pattern_begin_positions(size_t const read_len, uint64_t const pattern_size, uint64_t const overlap, functor_t && callback)
{
    assert(read_len >= pattern_size);
    assert(pattern_size >= overlap);

    size_t last_begin{0u};
    for (size_t i = 0; i <= read_len - pattern_size; i = i + pattern_size - overlap)
    {
        callback(i);
        last_begin = i;
    }

    if (last_begin < read_len - pattern_size)
    {
        // last pattern might have a smaller overlap to make sure the end of the read is covered
        callback(read_len - pattern_size);
    }
}

//-----------------------------
//
// Half open interval [begin; end) on the minimiser vector that shows minimisers belonging to the current pattern
// & the threshold for a local match.
//-----------------------------
struct pattern_bounds
{
    size_t begin_position;
    size_t end_position;
    size_t threshold;
};

//-----------------------------
//
// Match a pattern to the corresponding range of minimisers in the minimiser vector &
// find the threshold based on the number of minimisers in that range.
// NB: for winnowing minimisers patterns of the same length can contain a varying number of minimisers.
//
//-----------------------------
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

    size_t const minimiser_count = pattern.end_position - pattern.begin_position;

    pattern.threshold = thresholder.get(minimiser_count);

    return pattern;
}

//-----------------------------
//
// Count matching k-mers and return matching bins for one pattern.
//
//-----------------------------
template <typename binning_bitvector_t>
std::set<size_t> find_pattern_bins(pattern_bounds const & pattern, size_t const & bin_count, binning_bitvector_t const & counting_table)
{
    std::set<size_t> pattern_hits{};

    // counting vector for the current pattern
    seqan3::counting_vector<uint8_t> total_counts(bin_count, 0);

    for (size_t i = pattern.begin_position; i < pattern.end_position; i++)
        total_counts += counting_table[i];
    for (size_t current_bin = 0; current_bin < total_counts.size(); current_bin++)
    {
        auto &&count = total_counts[current_bin];
        if (count >= pattern.threshold)
        {
            // the result_set is a union of results from all patterns of a read
            pattern_hits.insert(current_bin);
        }
    }

    return pattern_hits;
}

struct local_prefilter_fn
{
    //-----------------------------
    //
    // Search a batch of queries in the IBF:
    // - the IBF is searched for local matches of length pattern size
    // - not all local match positions are considered
    // - overlap shows how much consequtive patterns overlap
    //
    //-----------------------------
    template <seqan3::data_layout ibf_data_layout>
    std::pair<std::vector<query_result>, std::map<size_t, std::vector<std::pair<std::string, std::vector<seqan3::dna4>>>>>
    operator()(
        std::span<query_record const> const & records,
        seqan3::interleaved_bloom_filter<ibf_data_layout> const & ibf,
        search_arguments const & arguments,
        raptor::threshold::threshold const & thresholder) const
    {
        // concurrent invocations of the membership agent are not thread safe
        // agent has to be created for each thread
        auto agent = ibf.membership_agent();
        size_t const bin_count = ibf.bin_count();

        std::vector<query_result> thread_query_result{}; // set of query results processed by one thread
        using query_list = std::vector<std::pair<std::string, std::vector<seqan3::dna4>>>;
        std::map<size_t, query_list> thread_bin_result{};

        // vector holding all the minimisers and their starting position for the read
        std::vector<std::tuple<uint64_t, size_t>> minimiser;

        auto minimiser_hash_adaptor = seqan3::views::minimiser_hash(
            arguments.shape,
            seqan3::window_size{arguments.window_size},
            seqan3::seed{adjust_seed(arguments.shape_weight)});

        for (query_record const & record : records)
        {
            std::string const & id = record.sequence_id;
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
            //	rows: each minimiser of read
            // 	columns: each bin of IBF
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

            std::set<size_t> sequence_hits{};
            pattern_begin_positions(seq.size(), arguments.pattern_size, arguments.overlap, [&](size_t const begin)
            {
                pattern_bounds const pattern = make_pattern_bounds(begin, arguments, window_span_begin, thresholder);
                std::set<size_t> const pattern_hits = find_pattern_bins(pattern, bin_count, counting_table);
                sequence_hits.insert(pattern_hits.begin(), pattern_hits.end());
            });

            for (size_t bin : sequence_hits)
            {
                thread_bin_result[bin].emplace_back(id, seq);
            }
            thread_query_result.emplace_back(id, std::move(sequence_hits));
        }


        return std::make_pair(std::move(thread_query_result), std::move(thread_bin_result));
    }
};

static constexpr local_prefilter_fn local_prefilter{};

} // namespace valik
