
#pragma once

#include <type_traits>

#include <dream_stellar/stellar_extension.hpp>
#include <dream_stellar/stellar_types.hpp>
#include <dream_stellar/utils/stellar_kernel_runtime.hpp>

namespace dream_stellar
{

///////////////////////////////////////////////////////////////////////////////
// Appends a segment of only error positions from align to queue.
template<typename TAlign, typename TPos, typename TScoreValue>
inline void
_appendNegativeSegment(TAlign const & align,
                       TPos & pos, TPos const len,
                       Score<TScoreValue> const & scoreMatrix,
                       String<Triple<TPos, TPos, TScoreValue> > & queue) {
    typedef Triple<TPos, TPos, TScoreValue> TMerger;
    TPos beginPos = pos;

    TScoreValue score = 0;
    while (pos < len) {
        if (isGap(row(align, 0), pos) || isGap(row(align, 1), pos)) {
            score += scoreGap(scoreMatrix);
        } else if (value(row(align, 0), pos) != value(row(align, 1), pos)) {
            score += scoreMismatch(scoreMatrix);
        } else {
            break;
        }
        ++pos;
    }
    if (pos == len) appendValue(queue, TMerger(beginPos, pos, MinValue<TScoreValue>::VALUE + 1));
    else appendValue(queue, TMerger(beginPos, pos, score));
}

///////////////////////////////////////////////////////////////////////////////
// Appends a segment of only matching positions from align to queue.
template<typename TAlign, typename TPos, typename TScoreValue>
inline void
_appendPositiveSegment(TAlign const & align,
                       TPos & pos, TPos const len,
                       Score<TScoreValue> const & scoreMatrix,
                       String<Triple<TPos, TPos, TScoreValue> > & queue) {
    if (pos == len) return;
    typedef Triple<TPos, TPos, TScoreValue> TMerger;
    TPos beginPos = pos;

    TScoreValue score = 0;
    while ((pos < len) &&
           (!isGap(row(align, 0), pos) &&
            !isGap(row(align, 1), pos) &&
            (value(row(align, 0), pos) == value(row(align, 1), pos)))) {
        score += scoreMatch(scoreMatrix);
        ++pos;
    }
    appendValue(queue, TMerger(beginPos, pos, score));
}

///////////////////////////////////////////////////////////////////////////////
// See Lemma 5 in Zhang et al., 1999.
template<typename TMerger>
inline bool
_negativeMerge(String<TMerger> & queue) {
    typedef typename Value<TMerger, 1>::Type TPos;
    TPos len = length(queue);
    if (len < 3) return false;

    TMerger cd = value(queue, len-1);
    TMerger bc = value(queue, len-2);
    TMerger ab = value(queue, len-3);

    if ((bc.i3 < 0) || (bc.i3 >= abs(_max(ab.i3, cd.i3)))) {
        return false;
    } else {
        String<TMerger> newMerger;
        appendValue(newMerger, TMerger(ab.i1, cd.i2, ab.i3 + bc.i3 + cd.i3));
        replace(queue, len-3, len, newMerger);

        return true;
    }
}

///////////////////////////////////////////////////////////////////////////////
// See Lemma 6 in Zhang et al., 1999.
template<typename TMerger>
inline bool
_positiveMerge(String<TMerger> & queue) {
    typedef typename Value<TMerger, 1>::Type TPos;
    TPos len = length(queue);
    if (len < 5) return false;

    TMerger ef = value(queue, len-1);
    TMerger de = value(queue, len-2);
    TMerger cd = value(queue, len-3);
    TMerger bc = value(queue, len-4);
    TMerger ab = value(queue, len-5);

    if ((cd.i3 >= 0) || (cd.i3 < _max(ab.i3, ef.i3))) {
        return false;
    } else {
        String<TMerger> newMerger;
        appendValue(newMerger, TMerger(bc.i1, de.i2, bc.i3 + cd.i3 + de.i3));
        replace(queue, len-4, len-1, newMerger);

        return true;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Implements the algorithm from Zhang et al. in Bioinformatics, 1999: "Post-processing long pairwise alignments".
// Splits an alignment into sub-alignments that contain no x-Drop.
template<typename TAlign, typename TScoreValue, typename TScoreValue1, typename TScoreValue2>
void
_splitAtXDrops(TAlign const & align,
               Score<TScoreValue> const & scoreMatrix,
               TScoreValue1 const scoreDropOff,
               TScoreValue2 const minScore,
               String<TAlign> & alignmentString) {
    typedef typename Position<Row<TAlign> >::Type TPos;
    typedef Triple<TPos, TPos, TScoreValue> TMerger;

    // initialization
    String<TMerger> queue;
    TPos pos = _min(toViewPosition(row(align, 0), beginPosition(row(align, 0))),
                    toViewPosition(row(align, 1), beginPosition(row(align, 1))));
    appendValue(queue, TMerger(pos, pos, MinValue<TScoreValue1>::VALUE + 1));

    TPos aliLength = _max(toViewPosition(row(align, 0), endPosition(row(align, 0))),
                          toViewPosition(row(align, 1), endPosition(row(align, 1))));
    TPos len;
    while ((pos < aliLength) || (length(queue) > 1)) {
        // construct useful tree
        if (!_negativeMerge(queue)) {
            if (!_positiveMerge(queue)) {
                _appendPositiveSegment(align, pos, aliLength, scoreMatrix, queue);
                _appendNegativeSegment(align, pos, aliLength, scoreMatrix, queue);
            }
        }

        // check for x-Drop
        len = length(queue);
        if ((len == 3) && (value(queue, 2).i3 < scoreDropOff * (-1))) {
            if (value(queue, 1).i3 >= minScore) {
                // create new sub-alignment
                TAlign ali(align);
                // std::cerr << "SEQ0 " << source(row(ali, 0)) << "\n";
                // std::cerr << "SEQ1 " << source(row(ali, 1)) << "\n";
                // std::cerr << "ROW0\n" << row(ali, 0) << "\nqueue[1].i1 == " << queue[1].i1 << ", queue[1].i2 == " << queue[1].i2 << "\n";
                // std::cerr << "ROW1\n" << row(ali, 1) << "\nqueue[1].i1 == " << queue[1].i1 << ", queue[1].i2 == " << queue[1].i2 <<
                        // "\n";
                setClippedBeginPosition(row(ali, 0), queue[1].i1 + clippedBeginPosition(row(align, 0)));
                setClippedBeginPosition(row(ali, 1), queue[1].i1 + clippedBeginPosition(row(align, 1)));
                setClippedEndPosition(row(ali, 0), queue[1].i2 + clippedBeginPosition(row(align, 0)));
                setClippedEndPosition(row(ali, 1), queue[1].i2 + clippedBeginPosition(row(align, 1)));

                // TPos begin0 = toSourcePosition(row(ali, 0), value(queue, 1).i1);
                // TPos begin1 = toSourcePosition(row(ali, 1), value(queue, 1).i1);
                // TPos end0 = toSourcePosition(row(ali, 0), value(queue, 1).i2);
                // TPos end1 = toSourcePosition(row(ali, 1), value(queue, 1).i2);
                // setClippedBeginPosition(row(ali, 0), begin0);
                // setClippedBeginPosition(row(ali, 1), begin1);
                // setBeginPosition(row(ali, 0), 0);
                // setBeginPosition(row(ali, 1), 0);
                // setClippedEndPosition(row(ali, 0), end0);
                // setClippedEndPosition(row(ali, 1), end1);

                // append sub-alignment
                appendValue(alignmentString, ali);
            }
            replace(queue, 0, 2, String<TMerger>());
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Conducts banded local alignment on swift hit (= computes eps-cores),
//  splits eps-cores at X-drops, and calls _extendAndExtract for extension of eps-cores
template<typename TSequence, typename TEpsilon, typename TSize, typename TDelta, typename TDrop,
         typename TOnAlignmentResultFn, bool bestLocalMethod>
void
allOrBestLocal(Segment<Segment<TSequence const, InfixSegment>, InfixSegment> const & infH,
               Segment<Segment<TSequence const, InfixSegment>, InfixSegment> const & infV,
               TEpsilon const eps,
               TSize const minLength,
               TDrop const xDrop,
               TDelta const delta,
               TOnAlignmentResultFn && onAlignmentResult,
               stellar_verification_time & verification_runtime,
               std::integral_constant<bool, bestLocalMethod>) {
    using TInfix = Segment<TSequence const, InfixSegment>;
    typedef Segment<TInfix, InfixSegment> TSegment;
    typedef typename StellarMatch<TSequence const, seqan2::CharString>::TAlign TAlign;

    TSize maxLength = 1000000000;
    if ((TSize)length(infH) > maxLength) {
        std::cerr << "Warning: SWIFT hit <" << beginPosition(infH) << "," << endPosition(infH);
        std::cerr << "> , <" << beginPosition(infV) << "," << endPosition(infV);
        std::cerr << "> too long. Verification skipped.\n" << std::flush;
        return;
    }

    // define a scoring scheme
    typedef int TScore;
    TScore match = 1;
    // large negative scoring scheme values lead to excessive seed extension
    TScore scoringSchemeLowerBound = -(TScore)1000;
    TScore mismatchIndel = scoringSchemeLowerBound;
    if (eps > -1/scoringSchemeLowerBound)   // avoid division by 0
        mismatchIndel = (TScore)_max((TScore) std::ceil(-1/eps) + 1, -(TScore)length(host(infH)));
    Score<TScore> scoreMatrix(match, mismatchIndel, mismatchIndel);
    TScore scoreDropOff = (TScore) _max((TScore) xDrop * (-mismatchIndel), MinValue<TScore>::VALUE + 1);

    // calculate minimal score for local alignments
    TEpsilon e = floor(eps*minLength);
    TSize minLength1 = _max(0, (TSize)ceil((e+1) / eps));
    TEpsilon e1 = floor(eps*minLength1);
    TScore minScore = _min((TSize)ceil((minLength-e) / (e+1)), (TSize)ceil((minLength1-e1) / (e1+1)));

    // diagonals for banded local alignment
    int64_t upperDiag = 0;
    int64_t lowerDiag = endPosition(infH) - (int64_t)endPosition(infV) - beginPosition(infH) + beginPosition(infV);

    std::cerr << "endPosition(infH)\t" << endPosition(infH) << '\n';
    std::cerr << "(int64_t)endPosition(infV)\t" << (int64_t)endPosition(infV) << '\n';
    std::cerr << "beginPosition(infH)\t" << beginPosition(infH) << '\n';
    std::cerr << "beginPosition(infV)\t" << beginPosition(infV) << '\n';
    std::cerr << "length infH\t" << endPosition(infH) - (int64_t)beginPosition(infH) << '\n';
    std::cerr << "infH\t" << infH << '\n';
    
    std::cerr << "length infV\t" << endPosition(infV) - (int64_t)beginPosition(infV) << '\n';
    std::cerr << "infV\t" << infV << '\n';

    std::cerr << "lowerDiag\t" << lowerDiag << '\n';

    if (beginPosition(infV) == 0) {
        if (endPosition(infV) == endPosition(host(infV))) {
            // TODO: is it possible to get a smaller band in this case?
            upperDiag = delta;
            lowerDiag = -(int64_t)delta;
        } else
            upperDiag = lowerDiag + delta;
    } else if (endPosition(infV) == endPosition(host(infV))) {
        lowerDiag = -(int64_t)delta;
    } 
    /*
    else if ((endPosition(infH) - (int64_t)beginPosition(infH)) > (endPosition(infV) - beginPosition(infV))) {
        upperDiag = delta;
        lowerDiag = -(int64_t)delta;
    }
    */
    
    std::cerr << "delta\t" << delta << '\n';
    std::cerr << "upperDiag\t" << upperDiag << '\n';
    std::cerr << "lowerDiag\t" << lowerDiag << '\n';

    // banded local alignment
    LocalAlignmentEnumerator<Score<TScore>, Banded> enumerator(scoreMatrix, lowerDiag, upperDiag, minScore);
    Align<TSegment> localAlign;
    resize(rows(localAlign), 2);
    assignSource(row(localAlign, 0), infH);
    assignSource(row(localAlign, 1), infV);

    while (true) {
        bool const has_next = verification_runtime.next_local_alignment_time.measure_time([&]()
        {
            return nextLocalAlignment(localAlign, enumerator);
        });

        if (!has_next)
            break;

    // while (localAlignment(localAlign, finder, scoreMatrix, minScore, lowerDiag, upperDiag, BandedWatermanEggert())) {

        // std::cerr << "localAlign == \n" << localAlign << "\n";

        // split local alignments containing an X-drop
        String<Align<TSegment> > seedAlignments;
        verification_runtime.split_at_x_drops_time.measure_time([&]()
        {
            _splitAtXDrops(localAlign, scoreMatrix, scoreDropOff, minScore, seedAlignments);
        });

        typename Iterator<String<Align<TSegment> > >::Type aliIt = begin(seedAlignments);
        while (aliIt != end(seedAlignments)) {
            // std::cerr << "aliIt == \n" << row(*aliIt, 0) << "\n" << row(*aliIt, 1) << "\n";
            // create alignment object for the complete sequences
            TAlign align;
            resize(rows(align), 2);
            setSource(row(align, 0), host(host(infH)));
            setSource(row(align, 1), host(host(infV)));

            // determine extension direction
            ExtensionDirection direction;
            if (length(seedAlignments) == 1) direction = EXTEND_BOTH;
            else if (aliIt == begin(seedAlignments)) direction = EXTEND_RIGHT;
            else if (aliIt == end(seedAlignments)-1) direction = EXTEND_LEFT;
            else direction = EXTEND_NONE;

            bool const extension_succeeded = verification_runtime.extension_time.measure_time([&]()
            {
                return _extendAndExtract(*aliIt, scoreDropOff, scoreMatrix, infH, infV, direction, minLength, eps, align, verification_runtime.extension_time);
            });

            // extend alignment and obtain longest contained eps-match
            if (!extension_succeeded) {
                aliIt++;
                continue;
            }

            // insert eps-match in matches string
            bool success = onAlignmentResult(align);
            if (!success)
                return;

            ++aliIt;
        }
        if (bestLocalMethod) break;
    }
}

} // namespace dream_stellar
