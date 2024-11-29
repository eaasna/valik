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

#ifndef SEQAN_HEADER_STELLAR_H
#define SEQAN_HEADER_STELLAR_H

#include <iostream>
#include <seqan/seeds.h>

#include <dream_stellar/stellar_types.hpp>
#include <dream_stellar/stellar_extension.hpp>
#include <dream_stellar/stellar_database_segment.hpp>
#include <dream_stellar/stellar_query_segment.hpp>
#include <dream_stellar/stellar_query_segment.tpp>
#include <dream_stellar/stellar_index.hpp>
#include <dream_stellar/utils/stellar_kernel_runtime.hpp>
#include <dream_stellar/verification/all_local.hpp>
#include <dream_stellar/verification/best_local.hpp>
#include <dream_stellar/verification/swift_hit_verifier.hpp>

#include <dream_stellar/diagnostics/print.hpp>

namespace dream_stellar
{

using namespace seqan2;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Checks whether two matches overlap *in seq2* and
//  whether the non-overlaping parts are shorter than minLength.
template<typename TMatch, typename TSize>
bool
checkOverlap(TMatch const & matchA, TMatch const & matchB, TSize const minLength) {
    // check id and orienation
    if (matchA.id != matchB.id || matchA.orientation != matchB.orientation) return false;
    if (matchA.id == TMatch::INVALID_ID || matchB.id == TMatch::INVALID_ID) return false;

    // check overlap in seq2
    if (matchA.begin2 >= matchB.begin2) {
        if (matchA.end2 >= matchB.end2) {
            // check length of non-overlapping parts of both matches
            if ((TSize)matchA.begin2 - (TSize)matchB.begin2 >= minLength &&
                (TSize)matchA.end2 - (TSize)matchB.end2 >= minLength) {
                return false;
            }
        }
        // check whether offset is the same in both sequences
        if (toViewPosition(matchA.row2, matchA.begin2) - toViewPosition(matchB.row2, matchB.begin2) !=
            toViewPosition(matchA.row1, matchA.begin1) - toViewPosition(matchB.row1, matchB.begin1)) {
            return false;
        }
    } else {
        if (matchA.end2 < matchB.end2) {
            // check length of non-overlapping parts of both matches
            if ((TSize)matchB.begin2 - (TSize)matchA.begin2 >= minLength &&
                (TSize)matchB.end2 - (TSize)matchA.end2 >= minLength) {
                return false;
            }
        }
        // check whether offset is the same in both sequences
        if (toViewPosition(matchB.row2, matchB.begin2) - toViewPosition(matchA.row2, matchA.begin2) !=
            toViewPosition(matchB.row1, matchB.begin1) - toViewPosition(matchA.row1, matchA.begin1)) {
            return false;
        }
    }
    return true;
}

template<typename TRow, typename TPosition>
TPosition
projectedPosition(TRow const & rowA, TRow const & rowB, TPosition pos)
{
    return toSourcePosition(rowB, toViewPosition(rowA, pos));
}

///////////////////////////////////////////////////////////////////////////////
// Checks all alignment columns of two overlapping matches.
// It is assumed that matchA.begin1 < matchB.begin1.
template<typename TMatch, typename TSize>
bool
_checkAlignColOverlap(TMatch const & matchA, TMatch const & matchB, TSize const minLength)
{
    TSize equalCols = 0;
    TSize diffCols = 0;

    for  (typename TMatch::TPos pos = matchB.begin1; pos < _min(matchA.end1, matchB.end1); ++pos)
    {
        if (projectedPosition(matchA.row1, matchA.row2, pos) == projectedPosition(matchB.row1, matchB.row2, pos))
            ++equalCols;
        else
            ++diffCols;
    }

    if (diffCols >= minLength) return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Marks matches that overlap in both sequences with a longer match as invalid.
template<typename TSequence, typename TId, typename TSize>
void maskOverlaps(String<StellarMatch<TSequence const, TId> > & matches, TSize const minLength)
{
    typedef StellarMatch<TSequence const, TId>              TMatch;
    typedef typename TMatch::TPos                           TPos;
    typedef typename Iterator<String<TMatch>, Rooted>::Type TIter;
    typedef typename Iterator<String<TSize>, Rooted>::Type  TOverlapIter;

    // sort matches by begin position in row0
    sortMatches(matches, LessPos<TMatch>());

    // iterate all matches
    TIter it = begin(matches);

    // list of matches that potentially overlap with the current match in row0 and
    // start earlier (including matches that overlap but have a unique part of at
    // least minLength) sorted by descending end positions
    String<TSize> overlaps;

    for (; it != end(matches); ++it)
    {
        if ((*it).id == TMatch::INVALID_ID) continue;

        TPos insertPos = 0;

        // iterate potentially overlapping matches
        TOverlapIter overlapIt = begin(overlaps);
        for (; overlapIt != end(overlaps); ++overlapIt)
        {
            TMatch & o = matches[*overlapIt];

            // determine position for inserting *it into overlaps after checking
            if ((*it).end1 < o.end1) insertPos++;

            // check if matches overlap in row0 - if not, then break
            if (o.end1 <= (*it).begin1) break;

            // check if unique parts of the two matches in row0 are longer than minLength - if yes, then continue
            if ((*it).begin1 - o.begin1 >= (TPos)minLength &&
                (*it).end1 > o.end1 && (*it).end1 - o.end1 >= (TPos)minLength) continue;

            // check if matches overlap in row1 - if not, then continue
            if (!checkOverlap(*it, o, minLength)) continue;

            // check exact alignment columns for overlap
            if (!_checkAlignColOverlap(o, *it, minLength)) continue;

            // set shorter match invalid
            if (length(*it) > length(o))
                o.id = TMatch::INVALID_ID;
            else
                (*it).id = TMatch::INVALID_ID;
        }

        // remove all matches from overlaps that end earlier than current match begins
        resize(overlaps, position(overlapIt));

        if ((*it).id != TMatch::INVALID_ID)
            insertValue(overlaps, insertPos, position(it));
    }

}

///////////////////////////////////////////////////////////////////////////////
// Removes matches that are marked as invalid, and then keeps only the numMatches best matches.
template<typename TSequence, typename TId, typename TSize>
void
compactMatches(String<StellarMatch<TSequence const, TId> > & matches, TSize const numMatches) {
    typedef StellarMatch<TSequence const, TId>                  TMatch;
    typedef typename Iterator<String<TMatch>, Standard>::Type   TIterator;

    // sort matches by length (and validity)
    sortMatches(matches, LessLength<TMatch>());

    // count valid matches
    TSize num = 0;

    TIterator it = begin(matches, Standard());
    TIterator itEnd = end(matches, Standard());

    for(; it != itEnd; ++it) {
        if ((*it).id != TMatch::INVALID_ID)
            ++num;
    }

    // keep only valid and longest matches
    resize(matches, _min(num, numMatches));
}

template<typename TMatch_>
inline bool
QueryMatches<TMatch_>::
removeOverlapsAndCompactMatches(size_t const disableThresh,
                                size_t const compactThresh,
                                size_t const minLength,
                                size_t const numMatches)
{
    if (this->disabled)
        return false;

    size_t const matchesCount = length(this->matches);

    if (matchesCount > disableThresh) {
        this->disabled = true;
        clear(this->matches);
        return false;
    }

    if (matchesCount <= compactThresh)
        return false;

    maskOverlaps(this->matches, minLength);      // remove overlaps and duplicates
    compactMatches(this->matches, numMatches);   // keep only the <numMatches> longest matches
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Appends a match to matches container and removes overlapping matches if threshold is reached.
template<typename TSource, typename TId, typename TSize, typename TSize1>
inline bool
_insertMatch(QueryMatches<StellarMatch<TSource const, TId> > & queryMatches,
             StellarMatch<TSource const, TId> const & match,
             TSize const minLength,
             TSize1 const disableThresh,
             TSize1 & compactThresh,
             TSize1 const numMatches) {

    appendValue(queryMatches.matches, match);

    // std::cerr << "Inserting match \n-------------\n" << match.row1 <<"\n" << match.row2 << "----------------\n";

    if (queryMatches.removeOverlapsAndCompactMatches(disableThresh, compactThresh, minLength, numMatches))
    {
        // raise compact threshold if many matches are kept
        if ((length(queryMatches.matches) << 1) > compactThresh)
            compactThresh += (compactThresh >> 1);
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Calls swift filter and verifies swift hits. = Computes eps-matches.
// A basic block for stellar
template<typename TAlphabet, typename TTag, typename TIsPatternDisabledFn, typename TOnAlignmentResultFn>
StellarComputeStatistics
_stellarKernel(StellarSwiftFinder<TAlphabet> & finder,  // iterate over database
               StellarSwiftPattern<TAlphabet> & pattern,    // holds the query and preprocessing info
               SwiftHitVerifier<TTag> & swiftVerifier,
               TIsPatternDisabledFn && isPatternDisabled,
               TOnAlignmentResultFn && onAlignmentResult,
               stellar_kernel_runtime & stellar_kernel_runtime) {
    StellarComputeStatistics statistics{};

    //!TODO: the pattern is non-const here. Where is it being updated?
    std::cout << "stellarKernel\n";
    auto const & index = host(pattern);
    std::cout << "BucketMap qgramCode length\t" << seqan2::length(index.bucketMap.qgramCode) << '\n'; 
    std::cout << "BucketMap qgramCode length\t" << seqan2::length(seqan2::indexDir(index)) << '\n';
    /*
    for (uint64_t i{1}; i < seqan2::length(index.bucketMap.qgramCode); i++)
    {
        uint64_t bucketBegin = index.bucketMap.qgramCode[i - 1];
        uint64_t bucketEnd = index.bucketMap.qgramCode[i];
        if (bucketBegin > bucketEnd)
        {
            std::cout << "Error: bucketBegin > bucketEnd\n" << bucketBegin << " > " << bucketEnd << "\n";
        }
    }
    */

    while (true) {

        bool const has_next = stellar_kernel_runtime.swift_filter_time.measure_time([&]()
        {
            return find(finder, pattern, swiftVerifier.eps_match_options.epsilon, swiftVerifier.eps_match_options.minLength);
        });

        if (!has_next)
            break;

        StellarDatabaseSegment<TAlphabet> databaseSegment
            = StellarDatabaseSegment<TAlphabet>::fromFinderMatch(infix(finder));

        ++statistics.numSwiftHits;
        statistics.totalLength += databaseSegment.size();
        statistics.maxLength = std::max<size_t>(statistics.maxLength, databaseSegment.size());

        if (isPatternDisabled(pattern)) continue;

        StellarQuerySegment<TAlphabet> querySegment
            = StellarQuerySegment<TAlphabet>::fromPatternMatch(pattern);

        ////Debug stuff:
        //std::cout << beginPosition(infix(finder)) << ",";
        //std::cout << endPosition(infix(finder)) << "  ";
        //std::cout << beginPosition(pattern) << ",";
        //std::cout << endPosition(pattern) << std::endl;

        // verification
        stellar_kernel_runtime.verification_time.measure_time([&]()
        {
            swiftVerifier.verify(
                databaseSegment,
                querySegment,
                pattern.bucketParams[0].delta + pattern.bucketParams[0].overlap,
                onAlignmentResult,
                stellar_kernel_runtime.verification_time);
        }); // measure_time
    }

    return statistics;
}

} // namespace dream_stellar

#endif
