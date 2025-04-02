
#pragma once

#include <seqan/seq_io.h>

#include <dream_stellar/stellar_sequence_segment.hpp>
#include <dream_stellar/stellar_types.hpp>

namespace dream_stellar
{

template <typename TAlphabet>
struct StellarDatabaseSegment : public StellarSequenceSegment<TAlphabet>
{
    using TBase = StellarSequenceSegment<TAlphabet>;

    using typename TBase::TInfixSegment;
    using TNestedFinderSegment = seqan2::Segment<TInfixSegment, seqan2::InfixSegment>;

    using TBase::TBase; // import constructor

    static StellarDatabaseSegment<TAlphabet> fromFinderMatch(TInfixSegment const & finderMatch)
    {
        seqan2::String<TAlphabet> const & underlyingDatabase = host(finderMatch);
        return {underlyingDatabase, seqan2::beginPosition(finderMatch), seqan2::endPosition(finderMatch)};
    }

    seqan2::String<TAlphabet> const & underlyingDatabase() const &
    {
        return this->underlyingSequence();
    }

    TNestedFinderSegment asFinderSegment() const
    {
        seqan2::String<TAlphabet> const & _database = underlyingDatabase();
        auto finderInfix = this->asInfixSegment();

        TInfixSegment const finderInfixSeq = infix(_database, 0, length(_database));
        TNestedFinderSegment finderSegment(finderInfixSeq,
            seqan2::beginPosition(finderInfix) - seqan2::beginPosition(_database),
            seqan2::endPosition(finderInfix) - seqan2::beginPosition(_database));
        assert(seqan2::length(finderInfix) == seqan2::length(finderSegment));
        return finderSegment;
    }
};

template <typename TAlphabet, typename TStorage>
TStorage _getDatabaseSegments(StringSet<String<TAlphabet>> & databases, StellarOptions const & options, bool const reverse = false)
{
    TStorage databaseSegments{};
    if (options.searchSegment)
    {
        if (length(databases[0]) < options.segmentEnd)
            throw std::runtime_error{"Segment end out of range"};

        if (options.segmentEnd <= options.segmentBegin)
            throw std::runtime_error{"Incorrect segment definition"};

        if (options.segmentEnd < options.minLength + options.segmentBegin)
            throw std::runtime_error{"Segment shorter than minimum match length"};

        if (reverse)
        {
            reverseComplement(databases[0]);
            databaseSegments.emplace_back(databases[0], length(databases[0]) - options.segmentEnd, length(databases[0]) - options.segmentBegin);
        }
        else
            databaseSegments.emplace_back(databases[0], options.segmentBegin, options.segmentEnd);
    }
    else
        for (auto & database : databases)
        {
            if (reverse)
                reverseComplement(database);

            if (length(database) >= options.minLength)
                databaseSegments.emplace_back(database, 0u, length(database));
        }

    return databaseSegments;
}

template <typename TAlphabet, typename TDatabaseSegment>
TDatabaseSegment _getDREAMDatabaseSegment(String<TAlphabet> const & sequenceOfInterest,
                                          StellarOptions const & options,
                                          bool const reverse = false)
{
    if (length(sequenceOfInterest) < options.segmentEnd)
        throw std::runtime_error{"Segment end out of range"};

    if (options.segmentEnd <= options.segmentBegin)
        throw std::runtime_error{"Incorrect segment definition"};

    if (options.segmentEnd < options.minLength + options.segmentBegin)
        throw std::runtime_error{"Segment shorter than minimum match length"};

    if (reverse)
        return TDatabaseSegment(sequenceOfInterest, length(sequenceOfInterest) - options.segmentEnd, length(sequenceOfInterest) - options.segmentBegin);

    return TDatabaseSegment(sequenceOfInterest, options.segmentBegin, options.segmentEnd);

}

} // namespace dream_stellar
