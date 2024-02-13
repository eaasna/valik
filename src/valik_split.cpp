#include <valik/split/split.hpp>

namespace valik::app
{

/**
 * @brief Function that divides reference or query database into partially overlapping segments.
 *
 * @param arguments Command line arguments.
 */
void valik_split(split_arguments & arguments)
{
    if (arguments.split_index)
        arguments.seg_count = adjust_bin_count(arguments.seg_count_in);
    else
    {
        //!TODO: read number of segments from reference metadata
        arguments.seg_count = arguments.seg_count_in;
    }

    metadata meta(arguments);
    meta.save(arguments.meta_out);

    if (!arguments.only_split)
    {
        // ==========================================
        // Parameter deduction
        // ==========================================
        auto space = param_space();
        std::vector<kmer_attributes> attr_vec;
        if (!read_fn_confs(attr_vec))
            precalc_fn_confs(attr_vec);

        search_pattern pattern(arguments.errors, arguments.pattern_size);
        if (arguments.split_index)
        {
            auto best_params = get_best_params(pattern, meta, attr_vec);
            arguments.kmer_size = best_params.k;
            kmer_attributes attr = attr_vec[arguments.kmer_size - std::get<0>(space.kmer_range)];

            if (arguments.verbose)
            {
                std::cout << "\n-----------Index build parameters-----------\n";
                std::cout << "db length " << meta.total_len << "bp\n";
                std::cout << "min local match length " << arguments.pattern_size << "bp\n";
                std::cout << "max error rate " << arguments.error_rate << "\n";
                std::cout << "kmer size " << std::to_string(arguments.kmer_size) << '\n';
    
                find_thresholds_for_kmer_size(meta, attr, arguments.error_rate);
            }
        }
        else
        {
            kmer_attributes attr = attr_vec[arguments.kmer_size - std::get<0>(space.kmer_range)];
            auto best_threshold = kmer_lemma_threshold(arguments.pattern_size, (size_t) arguments.kmer_size, (size_t) arguments.errors);
            if (arguments.verbose)
            {
                std::cout.precision(3);
                std::cout << "\n-----------Search parameters-----------\n";
                std::cout << "segment len " << std::to_string((uint64_t) std::round(meta.total_len / (double) meta.seg_count)) << '\n';
                std::cout << "min local match length " << arguments.pattern_size << "bp\n";
                std::cout << "error rate " << arguments.error_rate << "\n";
                std::cout << "kmer size " << std::to_string(arguments.kmer_size) << '\n';
            }
            
            metadata ref_meta = metadata(arguments.ref_meta_path);
            filtering_request request(pattern, ref_meta, meta);
            if (best_threshold <= 1)
            {
                best_threshold = request.find_heuristic_threshold(attr); 
                if (arguments.verbose)
                {
                    std::cout << "Apply heuristic threshold ";
                }
            }
            else
            {
                if (arguments.verbose)
                    std::cout << "Use k-mer lemma threshold ";
            }

            if (arguments.verbose)
            {
                param_set best_params(arguments.kmer_size, best_threshold, space);
                std::cout << best_threshold << '\n';
                std::cout << "FNR " <<  attr.fnr_for_param_set(pattern, best_params) << '\n';
                std::cout << "FPR " << request.fpr(best_params) << '\n';
            }
        }
    }

    if (arguments.write_out && !arguments.metagenome)
    {
        if (arguments.split_index)
            write_reference_segments(meta, arguments.db_file);
        else
            write_query_segments(meta, arguments.db_file);
    }
}

} // namespace valik::app
