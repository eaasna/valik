#include <utilities/threshold/find.hpp>

namespace valik
{

param_set get_best_params(search_pattern const & pattern,
                          metadata const & ref_meta,
                          std::vector<kmer_loss> const & attribute_vec, 
                          bool const verbose) 
{
    param_space space;
    param_set best_params(attribute_vec[0].k, 1, space);
    auto best_score = score(attribute_vec[0], pattern, best_params, ref_meta);

    std::vector<std::vector<double>> scores;
    std::vector<std::vector<double>> fn_rates;
    std::vector<std::vector<double>> fp_rates;
    scores.reserve(std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);
    fn_rates.reserve(std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);
    fp_rates.reserve(std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);

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

    if (verbose)
    {
        std::cout << "\n-----------FNR-----------\n";
        std::cout << '\t';
        for (size_t t{1}; t <= space.max_thresh; t++)
            std::cout << "t=" << t << '\t';
        std::cout << '\n';

        for (size_t i{0}; i < fn_rates.size(); i++)
        {
            std::cout << "k=" << std::get<0>(space.kmer_range) + i << '\t';
            for (auto & cell : fn_rates[i])
                std::cout << cell << '\t';
            std::cout << '\n';
        }

        std::cout << "\n-----------FP per pattern-----------\n";
        std::cout << '\t';
        for (size_t t{1}; t <= space.max_thresh; t++)
            std::cout << "t=" << t << '\t';
        std::cout << '\n';

        for (size_t i{0}; i < fp_rates.size(); i++)
        {
            std::cout << "k=" << std::get<0>(space.kmer_range) + i << '\t';
            for (auto & cell : fp_rates[i])
                std::cout << cell << '\t';
            std::cout << '\n';
        }
    
    }

    return best_params;
}

/**
 * @brief For a chosen kmer size and up to some maximum error rate find the best thresholds. 
 *        Consider different error rates up to max_err.
*/
search_kmer_profile find_thresholds_for_kmer_size(metadata const & ref_meta,
                                                  kmer_loss const attr, 
                                                  uint8_t const max_errors)
{
    param_space space;
    search_kmer_profile kmer_thresh{attr.k, ref_meta.pattern_size};
    for (uint8_t errors{0}; errors <= max_errors && errors <= space.max_errors; errors++)
    {
        search_pattern pattern(errors, ref_meta.pattern_size);
        search_kind search_type{LEMMA};
        auto best_params = param_set(attr.k, kmer_lemma_threshold(pattern.l, attr.k, errors), space);
        
        if ((best_params.t < 4) || 
            (1 - pow(1 - ref_meta.pattern_spurious_match_prob(best_params), 5e3) > 0.05))
        {
            search_type = HEURISTIC;
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
        
        uint64_t max_len = ref_meta.max_segment_len(best_params);
        if ((pattern.l * 100) > max_len)
        {
            search_type = STELLAR;
            kmer_thresh.add_error_rate(errors, {best_params, pattern, search_type});
        }
        else
        {
            double fnr = attr.fnr_for_param_set(pattern, best_params);
            double fp_per_pattern = ref_meta.pattern_spurious_match_prob(best_params);
            kmer_thresh.add_error_rate(errors, {best_params, pattern, search_type, fnr, fp_per_pattern, max_len});
        }
    }

    return kmer_thresh;
}

}   // namespace valik
