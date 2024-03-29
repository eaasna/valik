#pragma once

#include <string>
#include <vector>
#include <sstream>

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
};

} // namespace valik
