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
    else
    {
        // bin count is multiple of 64
        arguments.seg_count = adjust_bin_count(arguments.seg_count_in);
    }

    metadata meta(arguments);
    meta.save(arguments.meta_out);

    if (arguments.verbose)
    {
        std::cout << "\n-----------Preprocessing reference database-----------\n";
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
        std::vector<kmer_loss> attr_vec;
        if (!read_fn_confs(attr_vec))
            precalc_fn_confs(attr_vec);

        search_pattern pattern(arguments.errors, arguments.pattern_size);
        if (arguments.verbose)
        {
            std::cout << "\n-----------Local match definition-----------\n";
            std::cout << "min length " << arguments.pattern_size << "bp\n";
            std::cout << "max error rate " << arguments.error_rate << '\n';
        }
        
        auto best_params = get_best_params(pattern, meta, attr_vec, arguments.verbose);
        arguments.kmer_size = best_params.k;
        kmer_loss attr = attr_vec[arguments.kmer_size - std::get<0>(space.kmer_range)];

        search_kmer_profile search_profile = find_thresholds_for_kmer_size(meta, attr, arguments.errors);
        if (arguments.verbose)
            search_profile.print();

        std::filesystem::path search_profile_file{arguments.meta_out};
        search_profile_file.replace_extension("arg");
        search_profile.save(search_profile_file);
    }

    if (arguments.write_out && !arguments.metagenome)
    {
        write_reference_segments(meta, arguments.db_file);
    }
}

} // namespace valik::app
