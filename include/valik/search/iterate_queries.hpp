#pragma once

#include <valik/search/producer_threads_parallel.hpp>
#include <valik/search/search_time_statistics.hpp>
#include <valik/split/metadata.hpp>

#include <stellar/utils/stellar_app_runtime.hpp>
#include <stellar/io/import_sequence.hpp>

#include <seqan/seq_io.h>

#include <seqan3/core/debug_stream.hpp>

namespace valik::app
{

/**
 * @brief Function that sends chunks of queries to the prefilter which then writes shopping carts onto disk.
 *
 * @param arguments Command line arguments.
 * @param ibf Interleaved Bloom Filter.
 * @param thresholder Threshold for number of shared k-mers.
 * @param queue Shopping cart queue for load balancing between prefiltering and Stellar search.
 */
template <typename ibf_t, typename cart_queue_t>
void iterate_distributed_queries(search_arguments const & arguments,
                                 ibf_t const & ibf,
                                 raptor::threshold::threshold const & thresholder,
                                 cart_queue_t & queue)
{
    using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq>;
    std::vector<query_record> query_records{};
    seqan3::sequence_file_input<dna4_traits, fields> fin{arguments.query_file};
    for (auto &&chunked_records : fin | seqan3::views::chunk((1ULL << 20) * 10))
    {
        query_records.clear();
        for (auto && fasta_record: chunked_records)
            query_records.emplace_back(std::move(fasta_record.id()), std::move(fasta_record.sequence()));

        prefilter_queries_parallel(ibf, arguments, query_records, thresholder, queue);
    }
}

/**
 * @brief Function that creates a query record from each query sequence and sends it to Stellar search.
 *
 * @param ref_seg_count Number of reference segments i.e the distribution granularity.
 * @param arguments Command line arguments.
 * @param queue Shopping cart queue for sending queries over to Stellar search.
 */
template <typename count_t>
void iterate_all_queries(count_t const ref_seg_count,
                         search_arguments const & arguments,
                         cart_queue<shared_query_record<seqan2::String<seqan2::Dna>>> & queue)
{
    using TSequence = seqan2::String<seqan2::Dna>;
    using TId = seqan2::CharString;
    std::vector<shared_query_record<TSequence>> query_records{};
    constexpr uint64_t chunk_size = (1ULL << 20) * 10;

    seqan2::SeqFileIn inSeqs;
    if (!open(inSeqs, arguments.query_file.c_str()))
    {
        throw std::runtime_error("Failed to open " + arguments.query_file.string() + " file.");
    }

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    TSequence seq;
    TId id;
    size_t seqCount{0};
    for (; !atEnd(inSeqs); ++seqCount)
    {
        readRecord(id, seq, inSeqs);
        idsUnique &= stellar::_checkUniqueId(uniqueIds, id);

        auto query_ptr = std::make_shared<TSequence>(seq);
        std::vector<seqan3::dna4> seg_vec{};
        for (auto & c : *query_ptr)
        {
            seqan3::dna4 nuc;
            nuc.assign_char(c);
            seg_vec.push_back(nuc);
        }

        query_records.emplace_back(seqan2::toCString(id), std::move(seg_vec), seqan2::infix(*query_ptr, 0, seqan2::length(*query_ptr)), query_ptr);

        if (query_records.size() > chunk_size)
        {
            search_all_parallel<shared_query_record<seqan2::String<seqan2::Dna>>>(ref_seg_count, arguments, query_records, queue);
            query_records.clear();
        }
    }

    if (!idsUnique)
        std::cerr << "WARNING: Non-unique query ids. Output can be ambiguous.\n";

    search_all_parallel<shared_query_record<seqan2::String<seqan2::Dna>>>(ref_seg_count, arguments, query_records, queue);    
}

/**
 * @brief Function that creates short query records from fasta file input and sends them for prefiltering.
 *
 * @tparam ibf_t Interleaved Bloom Filter type.
 * @param arguments Command line arguments.
 * @param ibf Interleaved Bloom Filter of the reference database.
 * @param thresholder Threshold for number of shared k-mers.
 * @param queue Shopping cart queue for load balancing between Valik prefiltering and Stellar search.
 */
template <typename ibf_t>
void iterate_short_queries(search_arguments const & arguments,
                           ibf_t const & ibf,
                           raptor::threshold::threshold const & thresholder,
                           cart_queue<shared_query_record<seqan2::String<seqan2::Dna>>> & queue)
{
    using TSequence = seqan2::String<seqan2::Dna>;
    using TId = seqan2::CharString;
    std::vector<shared_query_record<TSequence>> query_records{};
    constexpr uint64_t chunk_size = (1ULL << 20) * 10;

    seqan2::SeqFileIn inSeqs;
    if (!open(inSeqs, arguments.query_file.c_str()))
    {
        throw std::runtime_error("Failed to open " + arguments.query_file.string() + " file.");
    }

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    TSequence seq;
    TId id;
    size_t seqCount{0};
    for (; !atEnd(inSeqs); ++seqCount)
    {
        readRecord(id, seq, inSeqs);
        idsUnique &= stellar::_checkUniqueId(uniqueIds, id);

        auto query_ptr = std::make_shared<TSequence>(seq);
        std::vector<seqan3::dna4> seg_vec{};
        for (auto & c : *query_ptr)
        {
            seqan3::dna4 nuc;
            nuc.assign_char(c);
            seg_vec.push_back(nuc);
        }

        query_records.emplace_back(seqan2::toCString(id), std::move(seg_vec), seqan2::infix(*query_ptr, 0, seqan2::length(*query_ptr)), query_ptr);

        if (query_records.size() > chunk_size)
        {
            prefilter_queries_parallel<shared_query_record<seqan2::String<seqan2::Dna>>>(ibf, arguments, query_records, thresholder, queue);
         
            query_records.clear();
        }
    }

    if (!idsUnique)
        std::cerr << "WARNING: Non-unique query ids. Output can be ambiguous.\n";

    prefilter_queries_parallel<shared_query_record<seqan2::String<seqan2::Dna>>>(ibf, arguments, query_records, thresholder, queue);
}

/**
 * @brief Function that creates split query records from fasta file input and sends them for prefiltering.
 *
 * @tparam ibf_t Interleaved Bloom Filter type.
 * @param arguments Command line arguments.
 * @param ibf Interleaved Bloom Filter of the reference database.
 * @param thresholder Threshold for number of shared k-mers.
 * @param queue Shopping cart queue for load balancing between Valik prefiltering and Stellar search.
 * @param meta Metadata table for split query segments.
 */
template <typename ibf_t>
void iterate_split_queries(search_arguments const & arguments,
                           ibf_t const & ibf,
                           raptor::threshold::threshold const & thresholder,
                           cart_queue<shared_query_record<seqan2::String<seqan2::Dna>>> & queue,
                           metadata & meta)
{
    using TSequence = seqan2::String<seqan2::Dna>;
    using TId = seqan2::CharString;
    std::vector<shared_query_record<TSequence>> query_records{};
    constexpr uint64_t chunk_size = (1ULL << 20) * 10;

    seqan2::SeqFileIn inSeqs;
    if (!open(inSeqs, arguments.query_file.c_str()))
    {
        throw std::runtime_error("Failed to open " + arguments.query_file.string() + " file.");
    }

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    TId id{};
    size_t seqCount{0};
    for (; !atEnd(inSeqs); ++seqCount)
    {
        TSequence seq{};
	    seqan3::debug_stream << "Read query record to split\n";
    	readRecord(id, seq, inSeqs);
        idsUnique &= stellar::_checkUniqueId(uniqueIds, id);

        auto query_ptr = std::make_shared<TSequence>(std::move(seq));

        for (auto const & seg : meta.segments_from_ind(seqCount))
        {
            seqan2::Segment<TSequence const, seqan2::InfixSegment> inf = seqan2::infixWithLength(*query_ptr, seg.start, seg.len);
            std::vector<seqan3::dna4> seg_vec{};
            for (auto & c : inf)
            {
                seqan3::dna4 nuc;
                nuc.assign_char(c);
                seg_vec.push_back(nuc);
            }

            query_records.emplace_back(seqan2::toCString(id), std::move(seg_vec), inf, query_ptr);

            if (query_records.size() > chunk_size)
            {
                prefilter_queries_parallel<shared_query_record<seqan2::String<seqan2::Dna>>>(ibf, arguments, query_records, thresholder, queue);
                query_records.clear();
            }
        }
    }

    if (!idsUnique)
        std::cerr << "WARNING: Non-unique query ids. Output can be ambiguous.\n";

    prefilter_queries_parallel<shared_query_record<seqan2::String<seqan2::Dna>>>(ibf, arguments, query_records, thresholder, queue);
}

}   // namespace valik::app
