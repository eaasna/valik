#include <utilities/threshold/find.hpp>

namespace valik
{

param_set get_best_params(search_pattern const & pattern,
                          metadata const & ref_meta,
                          std::vector<kmer_attributes> const & attribute_vec) 
{
    param_space space;
    param_set best_params(attribute_vec[0].k, 1, space);
    auto best_score = score(attribute_vec[0], pattern, best_params, ref_meta);

    std::vector<std::vector<double>> scores;
    std::vector<std::vector<double>> fn_rates;
    std::vector<std::vector<double>> fp_rates;
    scores.reserve(std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);
    fn_rates.reserve(std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);
    fp_rates.reserve(space.max_thresh);

    for (size_t i{0}; i < attribute_vec.size(); i++)
    {
        auto attr = attribute_vec[i];
        std::vector<double> kmer_scores;
        std::vector<double> kmer_fn;
        std::vector<double> kmer_fp;
        uint8_t k = attr.k;
        if (k > pattern.l)
            break;
        else
        {
            for (uint8_t t = 1; t <= space.max_thresh; t++)
            {
                param_set params(k, t, param_space());
                auto param_score = score(attr, pattern, params, ref_meta);
                kmer_scores.push_back(param_score);
                kmer_fn.push_back(attr.fnr_for_param_set(pattern, params));
                kmer_fp.push_back(ref_meta.pattern_spurious_match_prob(params));
                if (param_score <= best_score)
                {
                    best_score = param_score;
                    best_params = params;
                }
            }
            scores.push_back(kmer_scores);
            fn_rates.push_back(kmer_fn);
            fp_rates.push_back(kmer_fp);
        }
    }

    /*
    std::cout << '\t';
    for (size_t i{1}; i <= space.max_thresh; i++)
        std::cout << i << '\t';
    std::cout << '\n';

    for (size_t j{0}; j < fn_rates.size(); j++)
    {
        std::cout << std::get<0>(space.kmer_range) + j << '\t';
        for (auto & cell : fn_rates[j])
            std::cout << cell << '\t';
        std::cout << '\n';
    }

    std::cout << '\t';
    for (size_t i{1}; i <= space.max_thresh; i++)
        std::cout << i << '\t';
    std::cout << '\n';

    for (size_t j{0}; j < fp_rates.size(); j++)
    {
        std::cout << std::get<0>(space.kmer_range) + j << '\t';
        for (auto & cell : fp_rates[j])
            std::cout << cell << '\t';
        std::cout << '\n';
    }
    */

    return best_params;
}

/**
 * @brief For a chosen kmer size and up to some maximum error rate find the best thresholds. 
 *        Consider different error rates up to max_err.
*/
void find_thresholds_for_kmer_size(metadata const & ref_meta,
                                   kmer_attributes const attr, 
                                   double const max_err)
{
    param_space space;
    std::cout.precision(3);
    std::cout << "Recommended shared " << std::to_string(attr.k) << "-mer thresholds for different error rates\n";
    std::cout << "error_rate\tthreshold_kind\tthreshold\tFNR\tFP_per_pattern\tmax_query_seg_len\n";

    auto best_params = param_set(attr.k, space.max_thresh, space);
    for (uint8_t errors{1}; errors <= std::ceil(ref_meta.pattern_size * max_err) && errors <= space.max_errors; errors++)
    {
        search_pattern pattern(errors, ref_meta.pattern_size);
        std::cout << errors / (double) pattern.l << '\t';
        if (kmer_lemma_threshold(pattern.l, attr.k, errors) > 1)
        {
            best_params.t = kmer_lemma_threshold(pattern.l, attr.k, errors);
            std::cout << "kmer lemma\t";
        }
        else
        {
            std::cout << "heuristic\t";
            double best_score = pattern.l;
            for (uint8_t t{1}; t <= space.max_thresh; t++)
            {
                auto params = param_set(attr.k, t, space);
                auto param_score = score(attr, pattern, params, ref_meta);
                if (param_score <= best_score)
                {
                    best_params = params;
                    best_score = param_score;
                }
            }
        }

        std::cout << std::to_string(best_params.t) << '\t' << attr.fnr_for_param_set(pattern, best_params) << '\t'
                  << ref_meta.pattern_spurious_match_prob(best_params) << '\t' 
                  << std::to_string(ref_meta.max_segment_len(best_params)) << '\n';
    }
}

}   // namespace valik
