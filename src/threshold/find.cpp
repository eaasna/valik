#include <utilities/threshold/find.hpp>

namespace valik
{

/**
 * @brief For a pattern (with min length and max number of errors) find the best k-mer size and threshold for a given reference database.  
*/
param_set get_best_params(search_pattern const & pattern,
                          metadata const & ref_meta,
                          fn_confs const & fn_attr,
                          double const information_content, 
                          bool const verbose) 
{
    param_space space = fn_attr.space;

    if (pattern.l > space.max_len)
    {
        for (uint8_t k = space.max_k(); k >= space.min_k(); k--)
        {
            utilities::kmer kmer{k};
            if (kmer.lemma_threshold(pattern.l, pattern.e) > THRESH_LOWER)
                return param_set(kmer, kmer.lemma_threshold(pattern.l, pattern.e));
        }

        throw std::runtime_error{"Unable to deduce threshold for min_len=" + std::to_string(pattern.l) + " and errors=" + std::to_string(pattern.e)};
    }

    param_set best_params(space.min_k(), 1, space);
    auto best_score = score(fn_attr.get_kmer_loss(space.min_k()), pattern, best_params, ref_meta, PATTERNS_PER_SEGMENT, information_content);

    std::vector<std::vector<double>> scores;
    std::vector<std::vector<double>> fn_rates;
    std::vector<std::vector<double>> fp_rates;
    scores.reserve(space.max_k() - space.min_k() + 1);
    fn_rates.reserve(space.max_k() - space.min_k() + 1);
    fp_rates.reserve(space.max_k() - space.min_k() + 1);

    for (uint8_t k{space.min_k()}; k <= space.max_k(); k++)
    {
        const kmer_loss & attr = fn_attr.get_kmer_loss(k); 
        std::vector<double> kmer_scores;
        std::vector<double> kmer_fn;
        std::vector<double> kmer_fp;
        if (k > pattern.l)
            break;
        else
        {
            for (uint8_t t = 1; t <= space.max_thresh; t++)
            {
                param_set params(k, t, param_space());
                auto param_score = score(attr, pattern, params, ref_meta, PATTERNS_PER_SEGMENT, information_content);
                kmer_scores.push_back(param_score);
                kmer_fn.push_back(attr.fnr_for_param_set(pattern, params));
                kmer_fp.push_back(ref_meta.pattern_spurious_match_prob(params, information_content));
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
            std::cout << "k=" << space.min_k() + i << '\t';
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
            std::cout << "k=" << space.min_k() + i << '\t';
            for (auto & cell : fp_rates[i])
                std::cout << cell << '\t';
            std::cout << '\n';
        }
    
    }

    return best_params;
}

/**
 * @brief For a chosen kmer shape and up to some maximum error rate find the best thresholds. 
 *        Consider different error rates up to max_err.
*/
search_kmer_profile find_thresholds_for_kmer_size(metadata const & ref_meta,
                                                  kmer_loss attr, 
                                                  uint8_t const max_errors, 
                                                  double const information_content)
{
    constexpr param_space space{};
    search_kmer_profile kmer_thresh{attr.kmer, ref_meta.pattern_size};
    bool is_ungapped_shape = !attr.kmer.is_gapped();
    if (is_ungapped_shape)
    {
        for (uint8_t errors{0}; errors <= max_errors; errors++)
        {
            search_pattern pattern(errors, ref_meta.pattern_size);
            search_kind search_type{search_kind::LEMMA};

            if (attr.kmer.lemma_threshold(pattern.l, errors) > space.max_thresh)
            {
                auto best_params = param_set(attr.kmer, attr.kmer.lemma_threshold(pattern.l, errors));
                double fnr{0}; // == attr.fnr_for_param_set(pattern, best_params)
                double fp_per_pattern = ref_meta.pattern_spurious_match_prob(best_params, information_content);
                kmer_thresh.add_error_rate(errors, {best_params, pattern, search_type, fnr, fp_per_pattern});
            }
            else
            {
                auto best_params = param_set(attr.kmer, attr.kmer.lemma_threshold(pattern.l, errors), space);
                auto try_params = best_params;
                if ((best_params.t < THRESH_LOWER) ||  
                    segment_fpr(ref_meta.pattern_spurious_match_prob(best_params, information_content), PATTERNS_PER_SEGMENT) > FPR_UPPER)
                {
                    search_type = search_kind::STELLAR;
                    try_params.t++;
                    while ((try_params.t <= space.max_thresh) && 
                            (attr.fnr_for_param_set(pattern, try_params) < FNR_UPPER))
                    {
                        if (segment_fpr(ref_meta.pattern_spurious_match_prob(try_params, information_content), PATTERNS_PER_SEGMENT) <= FPR_UPPER)
                        {
                            best_params = std::move(try_params);
                            search_type = search_kind::HEURISTIC;
                        }
                        try_params.t++;
                    }
                }

                if (search_type == search_kind::STELLAR)
                {
                    kmer_thresh.add_error_rate(errors, {best_params, pattern, search_type});
                }
                else
                {
                    double fnr = attr.fnr_for_param_set(pattern, best_params);
                    double fp_per_pattern = ref_meta.pattern_spurious_match_prob(best_params, information_content);
                    kmer_thresh.add_error_rate(errors, {best_params, pattern, search_type, fnr, fp_per_pattern});
                }
            }
        }
    }
    else
    {
        for (uint8_t errors{0}; errors <= max_errors; errors++)
        {
            search_pattern pattern(errors, ref_meta.pattern_size);
            search_kind search_type{search_kind::HEURISTIC};

            if (attr.kmer.gapped_threshold(pattern.l, errors) > space.max_thresh)
            {
                auto best_params = param_set(attr.kmer, attr.kmer.gapped_threshold(pattern.l, errors));
                double fnr{0}; // == attr.fnr_for_param_set(pattern, best_params)
                double fp_per_pattern = ref_meta.pattern_spurious_match_prob(best_params, information_content);
                kmer_thresh.add_error_rate(errors, {best_params, pattern, search_type, fnr, fp_per_pattern});
            }
            else
            {
                auto best_params = param_set(attr.kmer, attr.kmer.gapped_threshold(pattern.l, errors), space);
                auto try_params = best_params;
                if ((best_params.t < THRESH_LOWER) ||  
                    segment_fpr(ref_meta.pattern_spurious_match_prob(best_params, information_content), PATTERNS_PER_SEGMENT) > FPR_UPPER)
                {
                    search_type = search_kind::STELLAR;
                    try_params.t++;
                    while ((try_params.t <= space.max_thresh) && 
                            (attr.fnr_for_param_set(pattern, try_params) < FNR_UPPER))
                    {
                        if (segment_fpr(ref_meta.pattern_spurious_match_prob(try_params, information_content), PATTERNS_PER_SEGMENT) <= FPR_UPPER)
                        {
                            best_params = std::move(try_params);
                            search_type = search_kind::HEURISTIC;
                        }
                        try_params.t++;
                    }
                }

                if (search_type == search_kind::STELLAR)
                {
                    kmer_thresh.add_error_rate(errors, {best_params, pattern, search_type});
                }
                else
                {
                    double fnr = attr.fnr_for_param_set(pattern, best_params);
                    double fp_per_pattern = ref_meta.pattern_spurious_match_prob(best_params, information_content);
                    kmer_thresh.add_error_rate(errors, {best_params, pattern, search_type, fnr, fp_per_pattern});
                }
            }
        }
    }

    return kmer_thresh;
}

}   // namespace valik
