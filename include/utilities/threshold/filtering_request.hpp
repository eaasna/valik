#pragma once

#include <utilities/threshold/find.hpp>
#include <utilities/threshold/search_pattern.hpp>

#include <valik/split/metadata.hpp>

namespace valik
{

/**
 * @brief The user requests filtering by setting the following parameters.
 * 
 * @param e Number of errors.
 * @param l Minimum length of local match i.e pattern length.
 * @param ref_meta Metadata for reference database.
 * @param query_meta Metadata for query database.
*/
struct filtering_request
{
    const search_pattern & pattern;
    const metadata & ref_meta;
    const metadata & query_meta;

    filtering_request(search_pattern const & pat, metadata const & ref, metadata const & query) : 
                      pattern(pat), ref_meta(ref), query_meta(query)
    {
        auto space = param_space();
        if (pat.e > space.max_errors)
            throw std::runtime_error{"error_count=" + std::to_string(pat.e) + " out of range [0, " + 
                                     std::to_string(space.max_errors) + "]"};
        if (pat.l > space.max_len)
            throw std::runtime_error{"min_len=" + std::to_string(pat.l) + " out of range [0, " + 
                                     std::to_string(space.max_len) + "]"};
    }

    /**
     * @brief An approximation of the probability of a query segment having a spurious match in a reference bin. 
    */
    double fpr(param_set const & params) const
    {
        double pattern_p = ref_meta.pattern_spurious_match_prob(params);
        uint64_t patterns_per_segment = std::round((query_meta.total_len / (double) query_meta.seg_count - pattern.l + 1) / (double) query_every);
        return segment_fpr(pattern_p, patterns_per_segment);
    }

};

}   // namespace valik
