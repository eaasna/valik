
#pragma once

#include <dream_stellar/stellar_sequence_segment.hpp>

namespace dream_stellar
{

template <typename TAlphabet>
struct StellarQuerySegment : public StellarSequenceSegment<TAlphabet>
{
    using TBase = StellarSequenceSegment<TAlphabet>;

    using typename TBase::TInfixSegment;
    using TNestedPatternSegment = seqan2::Segment<TInfixSegment, seqan2::InfixSegment>;

    using TBase::TBase; // import constructor

    template <typename TSwiftPattern>
    static StellarQuerySegment<TAlphabet> fromPatternMatch(TSwiftPattern const & swiftPattern);

    seqan2::String<TAlphabet> const & underlyingQuery() const &
    {
        return this->underlyingSequence();
    }

    TNestedPatternSegment asPatternSegment() const
    {
        seqan2::String<TAlphabet> const & _query = underlyingQuery();
        auto patternInfix = this->asInfixSegment();

        std::cerr << "asPatternSegment()" << '\n';
        /*std::cerr << "seqan2::beginPosition(patternInfix)\t" << seqan2::beginPosition(patternInfix) << '\n';
        std::cerr << "seqan2::beginPosition(patternInfixSeq)\t" << seqan2::beginPosition(patternInfixSeq) << '\n';
        std::cerr << "seqan2::beginPosition(patternInfix) - seqan2::beginPosition(patternInfixSeq)\t" << seqan2::beginPosition(patternInfix) - seqan2::beginPosition(patternInfixSeq) << '\n';
        std::cerr << "seqan2::endPosition(patternInfix)\t" << seqan2::endPosition(patternInfix) << '\n';
        std::cerr << "seqan2::endPosition(patternInfix) - seqan2::beginPosition(patternInfixSeq)\t" << seqan2::endPosition(patternInfix) - seqan2::beginPosition(patternInfixSeq) << '\n';
        */

        //using T = decltype(infix(_query, 0, length(_query)));
        //static_assert(std::same_as<typename T, void>);
        //static_assert(std::same_as<typename seqan2::Parameter_<T>::Type, void>);
        return {
            infix(_query, 0, length(_query)),
            seqan2::beginPosition(patternInfix),
            seqan2::endPosition(patternInfix)
        };
    }
};

} // namespace dream_stellar
