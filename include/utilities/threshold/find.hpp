#pragma once

#include <utilities/threshold/kmer_loss.hpp>
#include <utilities/threshold/search_kmer_profile.hpp>

#include <valik/split/metadata.hpp>

namespace valik
{

//!TODO: determine suitable parameters
constexpr double FPR_UPPER{0.1};
constexpr double FNR_UPPER{0.05};
constexpr uint8_t THRESH_LOWER{4};
constexpr size_t PATTERNS_PER_SEGMENT{5000};

/**
 * @brief The false positive probability of a query segment that contains partially overlapping patterns.
*/
inline double segment_fpr(double const pattern_p, size_t const patterns_per_segment)
{
    double none_match_p = pow(1 - pattern_p, patterns_per_segment);
    return std::min(1 - none_match_p, 1.0);
}

/**
* @brief The maximum length of a query segment that does not appear spuriously in reference bins. 
*/
inline uint64_t max_segment_len(double const pattern_p, size_t const pattern_size, uint8_t query_every)
{
    //!TODO: this const is the *minimum* number of patterns per segment
    // higher values are possible if pattern_p is small    
    return pattern_size + query_every * (PATTERNS_PER_SEGMENT - 1);
}

/**
* @brief Score of the objective function for a parameter set. Smaller values are better.
*/
inline double score(kmer_loss const & attr, 
                    search_pattern const & pattern, 
                    param_set const & params, 
                    metadata const & ref_meta,
                    size_t const PATTERNS_PER_SEGMENT)
{
    
    double none_match_p = pow(1 - ref_meta.pattern_spurious_match_prob(params), PATTERNS_PER_SEGMENT);
    double fpr = std::min(1 - none_match_p, 1.0);

    return attr.fnr_for_param_set(pattern, params) + fpr;
}

/**
 * @brief Look through the parameter space to find the best kmer size and threshold. 
*/
param_set get_best_params(search_pattern const & pattern, 
                          metadata const & ref_meta,
                          std::vector<kmer_loss> const & attribute_vec, 
                          bool const verbose);

/**
 * @brief For a chosen kmer size and some maximum error rate find the best threshold. 
*/
search_kmer_profile find_thresholds_for_kmer_size(metadata const & ref_meta,
                                             kmer_loss const attr, 
                                             uint8_t const max_errors);

}   // namespace valik
