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

    double information_content{0.35};
    if (!arguments.only_split)
    {
        // ==========================================
        // Parameter deduction
        // ==========================================
        auto space = param_space();
        fn_confs fn_attr = fn_confs(space);

        search_pattern pattern(arguments.errors, arguments.pattern_size);
        if (arguments.verbose)
        {
            std::cout << "\n-----------Local match definition-----------\n";
            std::cout << "min length " << arguments.pattern_size << "bp\n";
            std::cout << "max error rate " << arguments.error_rate << '\n';
        }
        
        if (arguments.kmer_size == std::numeric_limits<uint8_t>::max())
        {
            auto best_params = get_best_params(pattern, meta, fn_attr, information_content, arguments.verbose);
            arguments.kmer_size = best_params.kmer.size();
            arguments.shape = seqan3::shape{seqan3::ungapped(arguments.kmer_size)};
            arguments.shape_str = std::string(arguments.kmer_size, '1');
            arguments.shape_weight = arguments.kmer_size;
        }
        const kmer_loss & attr = fn_attr.get_kmer_loss(utilities::kmer{arguments.shape});

        search_kmer_profile search_profile = find_thresholds_for_kmer_size(meta, attr, arguments.errors, information_content);
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
