
#pragma once

// can we replace this include by some forwards?
// needed for seqan2::indexText, seqan2::needle
#include <seqan/index.h>

#include <dream_stellar/stellar_sequence_segment.hpp>

namespace dream_stellar
{

template <typename TAlphabet>
template <typename TSwiftPattern>
StellarQuerySegment<TAlphabet>
StellarQuerySegment<TAlphabet>::fromPatternMatch(TSwiftPattern const & swiftPattern)
{
    size_t const queryID = swiftPattern.curSeqNo;
    auto const & queryInfix = seqan2::getSequenceByNo(queryID, seqan2::indexText(seqan2::needle(swiftPattern)));
    static_assert(std::is_same_v<decltype(queryInfix), TInfixSegment const &>);
    auto const & underlyingQuery = host(queryInfix);
    static_assert(std::is_same_v<decltype(underlyingQuery), seqan2::String<TAlphabet> const &>);
    
    // Hit gets clipped in find_swift.hpp
    // needle(Pattern) only has access to the segment and not the underlying sequence length?
    // extend infix over underlying query sequence
    auto const queryInfixInfix = seqan2::infix(swiftPattern, queryInfix);

    int64_t hitBegin = swiftPattern.curBeginPos;
    int64_t hitEnd = swiftPattern.curEndPos;
    uint64_t hostLength = seqan2::length(underlyingQuery);

    //std::cerr << "hostLength\t" << hostLength << '\n';
    uint64_t segBegin = seqan2::beginPosition(queryInfixInfix);
    uint64_t segEnd = seqan2::endPosition(queryInfixInfix);

    assert(segEnd > segBegin);
    assert(hitEnd > 0);
    if (segEnd - segBegin < hitEnd)
        segEnd = std::clamp(segBegin + hitEnd, segEnd, hostLength);

    if (hitBegin < 0)
    {
        assert(hitBegin + segBegin > 0);
        segBegin = std::clamp(segBegin + hitBegin, 0UL, segBegin);
    }

    return {underlyingQuery, segBegin, segEnd};
}

} // namespace dream_stellar
