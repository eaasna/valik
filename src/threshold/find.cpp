#include <utilities/threshold/find.hpp>

namespace valik
{

/**
 * @brief Score of the objective function for a parameter set. Smaller values are better.
*/
double param_score(filtering_request const & request, param_set const & params, kmer_attributes const & attr)
{
    return attr.fnr_for_param_set(request, params) + request.fpr(params.k) / (double) params.t;
}

param_set get_best_params(param_space const & space, 
                          filtering_request const & request,
                          std::vector<kmer_attributes> const & attribute_vec) 
{
    param_set best_params(attribute_vec[0].k, 1, space);
    auto best_score = param_score(request, best_params, attribute_vec[0]);

    std::vector<std::vector<double>> scores;
    std::vector<std::vector<double>> fn_rates;
    std::vector<double> fp_rates;
    scores.reserve(std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);
    fn_rates.reserve(std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);
    fp_rates.reserve(std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);

    for (size_t i{0}; i < attribute_vec.size(); i++)
    {
        auto att = attribute_vec[i];
        std::vector<double> kmer_scores;
        std::vector<double> kmer_fn;
        size_t k = att.k;
        for (size_t t = 1; t <= space.max_thresh; t++)
        {
            param_set params(k, t, param_space());
            auto score = param_score(request, params, att);
            kmer_scores.push_back(score);
            kmer_fn.push_back(att.fnr_for_param_set(request, params));
            if (score < best_score)
            {
                best_score = score;
                best_params = params;
            }
        }
        scores.push_back(kmer_scores);
        fn_rates.push_back(kmer_fn);
        fp_rates.push_back(request.fpr(att.k));
    }

    return best_params;
}

/**
 * @brief For a chosen kmer size and error rate find the best threshold. 
*/
size_t find_threshold(param_space const & space, 
                      metadata const & meta,
                      search_arguments const & arguments,
                      kmer_attributes const att)
{

    filtering_request request(arguments.errors, arguments.pattern_size, meta.total_len, meta.seg_count);
    auto best_params = param_set(arguments.shape_size, space.max_thresh, space);
    double best_score = arguments.pattern_size;
    for (size_t t{1}; t <= space.max_thresh; t++)
    {
        auto params = param_set(att.k, t, space);
        auto score = param_score(request, params, att);
        if (score <= best_score)
        {
            best_params = params;
            best_score = score;
        }
    }

    std::cout.precision(3);
    if (arguments.verbose)
    {
        std::cout << "threshold " << best_params.t << '\n';
        std::cout << "FNR " <<  att.fnr_for_param_set(request, best_params) << '\n';
        std::cout << "FP_per_bin " << request.fpr(att.k) << '\n';
    }

    return best_params.t;
}

/**
 * @brief For a chosen kmer size and up to some maximum error rate find the best thresholds. 
*/
void find_thresholds_for_kmer_size(param_space const & space, 
                                  metadata const & meta,
                                  kmer_attributes const att, 
                                  double const max_err)
{
    std::cout.precision(3);
    std::cout << "Recommended shared " << att.k << "-mer thresholds for different error rates\n";
    std::cout << "error_rate\tthreshold_kind\tthreshold\tFNR\tFP_per_bin\n";

    auto best_params = param_set(att.k, space.max_thresh, space);
    for (size_t errors{1}; errors < (meta.segment_overlap() * max_err); errors++)
    {
        filtering_request request(errors, meta.segment_overlap(), meta.total_len, meta.seg_count);
        std::cout << errors / (double) request.l << '\t';
        if (att.fnr_for_param_set(request, best_params) == 0)
        {
            std::cout << "kmer lemma\t" << kmer_lemma_threshold(request.l, att.k, errors) << '\t' << 0;
        }
        else
        {
            std::cout << "heuristic\t";
            double best_score = request.l;
            for (size_t t{1}; t <= space.max_thresh; t++)
            {
                auto params = param_set(att.k, t, space);
                auto score = param_score(request, params, att);
                if (score <= best_score)
                {
                    best_params = params;
                    best_score = score;
                }
            }
            std::cout << best_params.t << '\t' << att.fnr_for_param_set(request, best_params);
        }

        std::cout << '\t' << request.fpr(att.k) << '\n';
    }
}


}   // namespace valik
