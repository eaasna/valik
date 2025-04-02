
#pragma once

#include <limits.h>

#include <seqan/index.h>
#include <dream_stellar/stellar_sequence_segment.hpp>

namespace dream_stellar
{

template <typename TAlphabet>
template <typename TSwiftPattern>
StellarQuerySegment<TAlphabet>
StellarQuerySegment<TAlphabet>::fromPatternMatch(TSwiftPattern const & swiftPattern)
{
    assert(swiftPattern.curEndPos > swiftPattern.curBeginPos); // SWIFT hit coordinates relative to query infix
    
    size_t const queryID = swiftPattern.curSeqNo;
    auto const & queryInfix = seqan2::getSequenceByNo(queryID, seqan2::indexText(seqan2::needle(swiftPattern)));
    static_assert(std::is_same_v<decltype(queryInfix), TInfixSegment const &>);
    auto const & underlyingQuery = host(queryInfix);
    static_assert(std::is_same_v<decltype(underlyingQuery), seqan2::String<TAlphabet> const &>);

    assert(seqan2::endPosition(queryInfix) > seqan2::beginPosition(queryInfix)); // Infix coordinates relative to query sequence
    assert(swiftPattern.curBeginPos < 0LL || seqan2::beginPosition(queryInfix) < (uint64_t) (std::numeric_limits<int64_t>::max() - swiftPattern.curBeginPos));
    assert((int64_t) seqan2::beginPosition(queryInfix) + swiftPattern.curBeginPos >= 0LL);  // swiftPattern.curBeginPos can be negative
    assert(seqan2::beginPosition(queryInfix) + swiftPattern.curEndPos <= seqan2::length(underlyingQuery)); 

    return {underlyingQuery, (size_t) seqan2::beginPosition(queryInfix) + swiftPattern.curBeginPos, seqan2::beginPosition(queryInfix) + swiftPattern.curEndPos};
}

} // namespace dream_stellar
