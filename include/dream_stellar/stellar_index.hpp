// ==========================================================================
//                    STELLAR - SwifT Exact LocaL AligneR
//                   http://www.seqan.de/projects/stellar/
// ==========================================================================
// Copyright (C) 2010-2012 by Birte Kehr
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STELLAR_INDEX_HPP
#define SEQAN_HEADER_STELLAR_INDEX_HPP

#include <seqan/index.h>

#include <span>

#include <dream_stellar/options/index_options.hpp>

namespace seqan2
{

using StellarSwiftLocal = Tag<struct StellarSwiftLocal_>;
template<>
struct Swift<StellarSwiftLocal>: public Swift<SwiftLocal>{};

/////////////////////////////////////////////////////////////
// Creates a new hit and appends it to the finders hit list
template <
typename THaystack,
typename TIndex,
typename TBucket,
typename TBucketParams,
typename TSize>
inline void _createHit(
    Finder<THaystack, Swift<StellarSwiftLocal> > & finder,
    Pattern<TIndex, Swift<StellarSwiftLocal> > & pattern,
    TBucket & bkt,
    TBucketParams & bucketParams,
    int64_t diag,
    TSize ndlSeqNo)
{
    typedef typename FindResult<Finder<THaystack, Swift<StellarSwiftLocal> >, Pattern<TIndex, Swift<StellarSwiftLocal> > >::Type THit;
    int64_t lastInc = (int64_t)bkt.lastIncrement - pattern.finderPosOffset;
    int64_t firstInc = (int64_t)bkt.firstIncrement - pattern.finderPosOffset;
    
    if(diag > lastInc)
    {
        // bucket is reused since last increment
        TSize reusePos = (bucketParams.reuseMask + 1) << bucketParams.logDelta;
        diag -= (int64_t)ceil((diag-lastInc)/(double)reusePos) * reusePos;
    }
    
    // determine width, height, and begin position in needle
    TSize width = lastInc - firstInc + length(pattern.shape);   // length of hit in haystack
    TSize height = width + bucketParams.delta + bucketParams.overlap;   // length of hit in needle
    
    int64_t ndlBegin = lastInc + length(pattern.shape) - diag - height;
    auto const & queryInfix = seqan2::getSequenceByNo(ndlSeqNo, seqan2::indexText(seqan2::needle(pattern)));
    
    int64_t queryInfixBegin = static_cast<int64_t>(seqan2::beginPosition(queryInfix));    
    // extend seed past the beginning of the infix segment
    if (ndlBegin + queryInfixBegin < 0LL)
    {
        ndlBegin = std::clamp<int64_t>(ndlBegin, -queryInfixBegin, 0);
        height = lastInc + length(pattern.shape) - diag - ndlBegin;
        assert(lastInc + length(pattern.shape) >= diag + ndlBegin); 
        // can not extend seed past the beginning of the host sequence
        if (height < width)
            return;
    }
    
    // extend seed past the end of the infix segment
    int64_t swiftHitEnd = queryInfixBegin + ndlBegin + height; 
    if (swiftHitEnd > 0 && (uint64_t) swiftHitEnd > seqan2::length(seqan2::host(queryInfix)))
    {
        height = std::clamp<TSize>(height, 0, std::max<int64_t>(0, (seqan2::length(seqan2::host(queryInfix)) - queryInfixBegin - ndlBegin)));
        // can not extend seed past the end of the host sequence
        if (height < width)
            return;
    }
    
    assert(height >= width); // the band width of alignment in verifySwiftHit is height - width
    // create the hit
    THit hit = {                //                              *
        firstInc,               // bucket begin in haystack     * *
        ndlSeqNo,               // needle seq. number           *   *
        ndlBegin,               // bucket begin in needle       *     *
        width,                  // bucket width (non-diagonal)    *   *
        height                  // bucket height                    * *
    };                          //                                    *
    
    // append the hit to the finders hit list
    appendValue(finder.hits, hit);
    }
}

namespace dream_stellar
{
using namespace seqan2;

template <typename TAlphabet, typename TString = String<TAlphabet>, typename TInfixSegment = seqan2::Segment<TString const, seqan2::InfixSegment>>
using StellarQGramStringSet = StringSet<TInfixSegment, Owner<> >;

template <typename TAlphabet>
using StellarQGramIndex = Index<StellarQGramStringSet<TAlphabet> const, IndexQGram<SimpleShape, OpenAddressing> >;

template <typename TAlphabet>
using StellarSwiftPattern = Pattern<StellarQGramIndex<TAlphabet>, Swift<StellarSwiftLocal> >;

template <typename TAlphabet>
using StellarSwiftFinder = Finder<Segment<String<TAlphabet> const, InfixSegment> const, Swift<StellarSwiftLocal> >;

template <typename TAlphabet>
struct StellarIndex
{
    using TSequence = seqan2::String<TAlphabet>;
    using TInfixSegment = seqan2::Segment<seqan2::String<TAlphabet> const, seqan2::InfixSegment>;
    using TQGramStringSet = StellarQGramStringSet<TAlphabet>;

    template <typename TSpec>
    StellarIndex(StringSet<TSequence, TSpec> const & queries, IndexOptions const & options)
        : StellarIndex{convertImplStringSet(queries), options}
    {}

    StellarIndex(StringSet<TInfixSegment> const & queries, IndexOptions const & options)
        : StellarIndex{convertSegmentStringSet(queries), options}
    {}

    StellarIndex(std::span<TInfixSegment> const & queries, IndexOptions const & options)
        : StellarIndex{convertImplSpan(queries), options}
    {}

    StellarIndex() = delete;
    StellarIndex(StellarIndex &&) = delete;
    StellarIndex(StellarIndex const &) = delete;
    StellarIndex & operator=(StellarIndex &&) = delete;
    StellarIndex & operator=(StellarIndex const &) = delete;

    void construct()
    {
        indexRequire(qgramIndex, QGramSADir());
    }

    StellarSwiftPattern<TAlphabet> createSwiftPattern()
    {
        return {qgramIndex};
    }

    static StellarQGramIndex<TAlphabet> & qgramIndexFromPattern(StellarSwiftPattern<TAlphabet> & pattern)
    {
        return host(pattern);
    }

    static TQGramStringSet const & sequencesFromPattern(StellarSwiftPattern<TAlphabet> & pattern)
    {
        return sequencesFromQGramIndex(qgramIndexFromPattern(pattern));
    }

    static TQGramStringSet const & sequencesFromQGramIndex(StellarQGramIndex<TAlphabet> & index)
    {
        return indexText(index);
    }

    TQGramStringSet const dependentQueries;
    StellarQGramIndex<TAlphabet> qgramIndex;

private:
    template <typename TOtherQGramStringSet, typename = std::enable_if_t<std::is_same_v<TOtherQGramStringSet, TQGramStringSet>>>
    StellarIndex(TOtherQGramStringSet && queries, IndexOptions const & options)
        : dependentQueries{std::forward<TOtherQGramStringSet>(queries)}, qgramIndex{dependentQueries}
    {
        resize(indexShape(qgramIndex), options.qGram);
        cargo(qgramIndex).abundanceCut = options.qgramAbundanceCut;
    }

    template <typename TSpec>
    static TQGramStringSet convertImplStringSet(StringSet<TSequence, TSpec> const & queries)
    {
        TQGramStringSet dependentQueries;
        for (TSequence const & query: queries)
            seqan2::appendValue(dependentQueries, seqan2::infix(query, 0, seqan2::length(query)));

        return dependentQueries;
    }

    static TQGramStringSet convertSegmentStringSet(StringSet<TInfixSegment> const & queries)
    {
        TQGramStringSet dependentQueries = queries;
        return dependentQueries;
    }

    static TQGramStringSet convertImplSpan(std::span<TInfixSegment> const & queries)
    {
        TQGramStringSet dependentQueries;
        for (TInfixSegment const & query: queries)
            seqan2::appendValue(dependentQueries, query);

        return dependentQueries;
    }
};

} // namespace dream_stellar


namespace seqan2 {

template <typename TAlphabet>
struct Cargo<::dream_stellar::StellarQGramIndex<TAlphabet>>
{
    typedef struct
    {
        double      abundanceCut;
    } Type;
};

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TAlphabet>
inline bool _qgramDisableBuckets(::dream_stellar::StellarQGramIndex<TAlphabet> & index)
{
    using TIndex = ::dream_stellar::StellarQGramIndex<TAlphabet>;
    using TDir = typename Fibre<TIndex, QGramDir>::Type;
    using TDirIterator = typename Iterator<TDir, Standard>::Type;
    using TSize = typename Value<TDir>::Type;

    TDir & dir   = indexDir(index);
    bool result  = false;
    unsigned counter = 0;
    TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
    if (thresh < 100)
        thresh = 100;

    TDirIterator it = begin(dir, Standard());
    TDirIterator itEnd = end(dir, Standard());
    for (; it != itEnd; ++it)
        if (*it > thresh)
        {
            *it = (TSize) - 1;
            result = true;
            ++counter;
        }

    if (counter > 0)
        std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

    return result;
}

} // namespace seqan2

#endif
