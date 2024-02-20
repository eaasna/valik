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
    if (arguments.only_split)
    {
        // use user parameter input
        arguments.seg_count = arguments.seg_count_in;
    }
    else if (arguments.split_index)
    {
        // bin count is multiple of 64
        arguments.seg_count = adjust_bin_count(arguments.seg_count_in);
    }
    else
    {
        std::filesystem::path kmer_thresh_file{arguments.ref_meta_path};
        kmer_thresh_file.replace_extension("arg");
        sharg::input_file_validator argument_input_validator{{"arg"}};
        argument_input_validator(kmer_thresh_file);
        kmer_thresholds kmer_thresh{kmer_thresh_file};
        arguments.max_segment_len= kmer_thresh.error_table[arguments.errors].max_segment_len;
        // seg_count is inferred in metagenome constructor
    }

    metadata meta(arguments);
    meta.save(arguments.meta_out);

    if (arguments.verbose)
    {
        if (arguments.split_index)
            std::cout << "\n-----------Preprocessing reference database-----------\n";
        else
            std::cout << "\n-----------Preprocessing queries-----------\n";

        std::cout << "database size " << meta.total_len << "bp\n";
        std::cout << "segment count " << meta.seg_count << '\n';
        std::cout << "segment len " << std::to_string((uint64_t) std::round(meta.total_len / (double) meta.seg_count)) << "bp\n";
    }

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
        if (arguments.verbose)
        {
            std::cout << "\n-----------Local match definition-----------\n";
            std::cout << "min length " << arguments.pattern_size << "bp\n";
            std::cout << "max error rate " << arguments.error_rate << '\n';

        }
        if (arguments.split_index)
        {
            auto best_params = get_best_params(pattern, meta, attr_vec);
            arguments.kmer_size = best_params.k;
            kmer_attributes attr = attr_vec[arguments.kmer_size - std::get<0>(space.kmer_range)];

            kmer_thresholds kmer_thresh = find_thresholds_for_kmer_size(meta, attr, arguments.error_rate);
            if (arguments.verbose)
                kmer_thresh.print();

            std::filesystem::path kmer_thresh_file{arguments.meta_out};
            kmer_thresh_file.replace_extension("arg");
            kmer_thresh.save(kmer_thresh_file);
        }
        else
        {
            //!TODO: Move query split into search
            std::filesystem::path kmer_thresh_file{arguments.ref_meta_path};
            kmer_thresh_file.replace_extension("arg");
            kmer_thresholds kmer_thresh{kmer_thresh_file};
            error_threshold error_thresh = kmer_thresh.get_error_thresh(arguments.errors);
            if (error_thresh.pattern != pattern)
                throw std::runtime_error("Incompatible search pattern");

            arguments.kmer_size = kmer_thresh.k;            
            kmer_attributes attr = attr_vec[arguments.kmer_size - std::get<0>(space.kmer_range)];

            metadata ref_meta = metadata(arguments.ref_meta_path);
            filtering_request request(pattern, ref_meta, meta);

            if (request.fpr(error_thresh.params) > 0.2)
                std::cerr << "WARNING: Prefiltering will be inefficient for a high error rate.\n";

            if (arguments.verbose)
            {
                std::cout.precision(3);
                std::cout << "\n-----------Search parameters-----------\n";
                std::cout << "kmer size " << std::to_string(error_thresh.params.k) << '\n';
                if (error_thresh.is_heuristic)
                    std::cout << "heuristic ";
                else
                    std::cout << "k-mer lemma ";
                std::cout << "threshold ";
                std::cout << std::to_string(error_thresh.params.t) << '\n';

                std::cout << "FNR " << attr.fnr_for_param_set(error_thresh.pattern, error_thresh.params) << '\n';
                std::cout << "FPR " << request.fpr(error_thresh.params) << '\n';
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
