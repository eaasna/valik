#pragma once

#include <string>
#include <vector>
#include <sstream>

#include <valik/split/metadata.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan/sequence.h>

namespace valik
{

/**
 * @brief Each record owns the underlying data.
 *
 */
struct query_record
{
    std::string sequence_id;
    std::vector<seqan3::dna4> sequence;
};

/**
 * @brief Query sequence resources are shared between records.
 *
 */
template <typename TSequence>
struct shared_query_record
{
    std::string sequence_id;
    std::vector<seqan3::dna4> sequence;
    seqan2::Segment<TSequence const, seqan2::InfixSegment> querySegment;
    std::shared_ptr<TSequence> underlyingData;

    shared_query_record(TSequence seq, std::string id) : sequence_id(std::move(id))
    {
        // make_shared returns a newly allocated object
        auto query_ptr = std::make_shared<TSequence>(std::move(seq));
        std::vector<seqan3::dna4> seg_vec{};
        for (auto & c : *query_ptr)
        {
            seqan3::dna4 nuc;
            nuc.assign_char(c);
            seg_vec.push_back(nuc);
        }

        sequence = std::move(seg_vec);
        querySegment = seqan2::infix(*query_ptr, 0, seqan2::length(*query_ptr));
        underlyingData = query_ptr;
    }

    shared_query_record(std::string id, metadata::segment_stats const & seg, std::shared_ptr<TSequence> query_ptr) : sequence_id(std::move(id))
    {
        seqan2::Segment<TSequence const, seqan2::InfixSegment> inf = seqan2::infixWithLength(*query_ptr, seg.start, seg.len);
        std::vector<seqan3::dna4> seg_vec{};
        for (auto & c : inf)
        {
            seqan3::dna4 nuc;
            nuc.assign_char(c);
            seg_vec.push_back(nuc);
        }

        sequence = std::move(seg_vec); 
        querySegment = inf;
        underlyingData = query_ptr;
    }
};

} // namespace valik
