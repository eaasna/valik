#pragma once

#include <utilities/threshold/kmer_attributes.hpp>
#include <utilities/threshold/kmer_thresholds.hpp>

#include <valik/split/metadata.hpp>

namespace valik
{

/**
* @brief Score of the objective function for a parameter set. Smaller values are better.
*/
inline double score(kmer_attributes const & attr, search_pattern const & pattern, param_set const & params, metadata const & ref_meta)
{
    return attr.fnr_for_param_set(pattern, params) + ref_meta.pattern_spurious_match_prob(params);
}

/**
 * @brief Look through the parameter space to find the best kmer size and threshold. 
*/
param_set get_best_params(search_pattern const & pattern, 
                          metadata const & ref_meta,
                          std::vector<kmer_attributes> const & attribute_vec);

/**
 * @brief For a chosen kmer size and some maximum error rate find the best threshold. 
*/
kmer_thresholds find_thresholds_for_kmer_size(metadata const & ref_meta,
                                             kmer_attributes const attr, 
                                             double const max_err);

}   // namespace valik
