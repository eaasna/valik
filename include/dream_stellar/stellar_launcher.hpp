#pragma once

#include <dream_stellar/query_id_map.hpp>
#include <dream_stellar/stellar.hpp>

namespace dream_stellar
{

/**
 * @return Any query has valid matches.
 */
template <typename TAlphabet, typename TId>
bool _postproccessQueryMatches(bool const databaseStrand, uint64_t const & refLen,
                               StellarOptions const & options,
                               StringSet<QueryMatches<StellarMatch<String<TAlphabet> const, TId> > > & matches,
                               std::vector<size_t> & disabledQueryIDs)
{
    using TSequence = String<TAlphabet>;

    bool foundMatches{false};
    for (size_t queryID = 0; queryID < length(matches); ++queryID)
    {
        QueryMatches<StellarMatch<TSequence const, TId>> & queryMatches = value(matches, queryID);

        queryMatches.removeOverlapsAndCompactMatches(options.disableThresh,
                                                     /*compactThresh*/ 0,
                                                     options.minLength,
                                                     options.numMatches);

        if (queryMatches.disabled)
            disabledQueryIDs.push_back(queryID);
        
        if (seqan2::length(queryMatches.matches) > 0)
            foundMatches = true;
    }

    if (databaseStrand)
        _postproccessLengthAdjustment(refLen, matches);
    
    return foundMatches;
}

template <typename TAlphabet, typename TId = CharString>
struct StellarLauncher
{
    template <typename visitor_fn_t>
    static constexpr StellarComputeStatistics _verificationMethodVisit(
        StellarVerificationMethod verificationMethod,
        visitor_fn_t && visitor_fn
    )
    {
        if (verificationMethod == StellarVerificationMethod{AllLocal{}})
            return visitor_fn(AllLocal());
        else if (verificationMethod == StellarVerificationMethod{BestLocal{}})
            return visitor_fn(BestLocal());
        return StellarComputeStatistics{};
    }

    static StellarComputeStatistics
    search_and_verify(
        StellarDatabaseSegment<TAlphabet> const databaseSegment,
        TId const & databaseID,
        bool const databaseStrand,
        StellarOptions & localOptions, // localOptions.compactThresh is out-param
        StellarSwiftPattern<TAlphabet> & localSwiftPattern,
        dream_stellar::stellar_kernel_runtime & strand_runtime,
        StringSet<QueryMatches<StellarMatch<String<TAlphabet> const, TId> > > & localMatches
    )
    {
        using TSequence = String<TAlphabet>;

        auto getQueryMatches = [&](auto const & pattern) -> QueryMatches<StellarMatch<TSequence const, TId> > &
        {
            return value(localMatches, pattern.curSeqNo);
        };

        auto isPatternDisabled = [&](StellarSwiftPattern<TAlphabet> & pattern) -> bool {
            QueryMatches<StellarMatch<TSequence const, TId> > & queryMatches = getQueryMatches(pattern);
            return queryMatches.disabled;
        };

        auto onAlignmentResult = [&](auto & alignment) -> bool {
            QueryMatches<StellarMatch<TSequence const, TId> > & queryMatches = getQueryMatches(localSwiftPattern);

            StellarMatch<TSequence const, TId> match(alignment, databaseID, databaseStrand);
            length(match);  // DEBUG: Contains assertion on clipping.

            // success
            return _insertMatch(
                queryMatches,
                match,
                localOptions.minLength,
                localOptions.disableThresh,
                // compactThresh is basically an output-parameter; will be updated in kernel and propagated back
                // outside of this function, the reason why StellarOptions can't be passed as const to this function.
                //!TODO: We might want to make this tied to the QueryMatches itself, as it should know then to consolidate
                // the matches
                localOptions.compactThresh,
                localOptions.numMatches
            );
        };

        // finder
        StellarSwiftFinder<TAlphabet> swiftFinder(databaseSegment.asInfixSegment(), localOptions.minRepeatLength, localOptions.maxRepeatPeriod);

        StellarComputeStatistics statistics = _verificationMethodVisit(
            localOptions.verificationMethod,
            [&]<typename TTag>(TTag) -> StellarComputeStatistics
            {
                SwiftHitVerifier<TTag> swiftVerifier
                {
                    STELLAR_DESIGNATED_INITIALIZER(.eps_match_options = , localOptions),
                    STELLAR_DESIGNATED_INITIALIZER(.verifier_options = , localOptions),
                };

                return _stellarKernel(swiftFinder, localSwiftPattern, swiftVerifier, isPatternDisabled, onAlignmentResult, strand_runtime);
            });

        return statistics;
    }
};

}   // namespace dream_stellar
