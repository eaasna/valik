#include <valik/argument_parsing/build.hpp>

namespace valik::app
{

void init_build_parser(sharg::parser & parser, build_arguments & arguments)
{
    init_shared_meta(parser);
    parser.add_option(arguments.out_path,
                      sharg::config{.short_id = '\0',
                      .long_id = "output",
                      .description = "Provide an output filepath.",
                      .required = true,
                      .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {}}});
    parser.add_option(arguments.ref_meta_path,
                    sharg::config{.short_id = '\0',
                    .long_id = "ref-meta",
                    .description = "Path to reference metadata created by valik split.",
                    .required = true,
                    .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.threads,
                    sharg::config{.short_id = '\0',
                    .long_id = "threads",
                    .description = "Choose the number of threads.",
                    .validator = positive_integer_validator{}});
    parser.add_flag(arguments.fast,
                      sharg::config{.short_id = '\0',
                      .long_id = "fast",
                      .description = "Build the index in fast mode when few false negatives can be tolerated in the following search."});
    parser.add_flag(arguments.manual_parameters,
                      sharg::config{.short_id = '\0',
                      .long_id = "without-parameter-tuning",
                      .description = "Do not read parameters from metadata.",
                      .advanced = true});
    /////////////////////////////////////////
    // Advanced options
    /////////////////////////////////////////
    parser.add_option(arguments.window_size,
                      sharg::config{.short_id = 'w',
                      .long_id = "window",
                      .description = "Choose the window size.",
                      .advanced = true,
                      .validator = positive_integer_validator{}});
    parser.add_option(arguments.kmer_size,
                      sharg::config{.short_id = 'k',
                      .long_id = "kmer",
                      .description = "Choose the kmer size.",
                      .advanced = true,
                      .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(arguments.hash,
                      sharg::config{.short_id = '\0',
                      .long_id = "hash",
                      .description = "Choose the number of hashes.",
                      .advanced = true,
                      .validator = sharg::arithmetic_range_validator{1, 5}});
    parser.add_option(arguments.size,
                      sharg::config{.short_id = '\0',
                      .long_id = "size",
                      .description = "Choose the size of the resulting IBF.",
                      .advanced = true,
                      .validator = size_validator{"\\d+\\s{0,1}[k,m,g,t,K,M,G,T]"}});
    parser.add_flag(arguments.verbose,
                    sharg::config{.short_id = '\0',
                    .long_id = "verbose",
                    .description = "Print verbose output.",
                    .advanced = true});    
    
    parser.add_subsection("Processing options");
    parser.add_option(arguments.kmer_count_min_cutoff,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer-count-min",
                                    .description = "Only store k-mers with at least (>=) x occurrences. "
                                                   "Mutually exclusive with --use-filesize-dependent-cutoff.",
                                    .validator = sharg::arithmetic_range_validator{0, 254}});
    parser.add_option(arguments.kmer_count_max_cutoff,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer-count-max",
                                    .description = "Only store k-mers with no more than (<=) x occurrences. "
                                                   "Mutually exclusive with --use-filesize-dependent-cutoff.",
                                    .validator = sharg::arithmetic_range_validator{1, 254}});
    parser.add_flag(arguments.use_filesize_dependent_cutoff,
                    sharg::config{.short_id = '\0',
                                  .long_id = "use-filesize-dependent-cutoff",
                                  .description = "Apply cutoffs from Mantis(Pandey et al., 2018). "
                                                 "Mutually exclusive with --kmer-count-cutoff."});
}

void run_build(sharg::parser & parser)
{
    build_arguments arguments{};
    init_build_parser(parser, arguments);
    try_parsing(parser);

    if (!arguments.manual_parameters)
    {
        std::filesystem::path search_profile_file{arguments.ref_meta_path};
        search_profile_file.replace_extension("arg");
        sharg::input_file_validator argument_input_validator{{"arg"}};
        argument_input_validator(search_profile_file);
        search_kmer_profile search_profile{search_profile_file};

        if (parser.is_option_set("kmer"))
        {
            seqan3::debug_stream << "WARNING: kmer size k=" << arguments.kmer_size << " will be updated to " << search_profile.k 
                                 << ". Set --without-parameter-tuning to force manual input.";
        }
        if (parser.is_option_set("window"))
        {
            seqan3::debug_stream << "WARNING: window size w=" << arguments.window_size << " will be updated. " 
                                 << "Set --without-parameter-tuning to force manual input.";
        }
        
        arguments.kmer_size = search_profile.k;
        arguments.window_size = search_profile.k;
        arguments.shape = search_profile.shape;
    }
    else
    {
        arguments.shape = seqan3::shape(seqan3::ungapped(arguments.kmer_size));
    }
    arguments.shape_weight = arguments.shape.count();

    if (arguments.kmer_size > arguments.window_size)
    {
        throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};
    }

    // ==========================================
    // Process bin_path:
    // if building from clustered sequences each line in input corresponds to a bin
    // if building from overlapping segments all sequences in one reference file
    // ==========================================
    auto sequence_file_validator{bin_validator{}.sequence_file_validator};

    metadata meta(arguments.ref_meta_path);
    arguments.bins = meta.seg_count;
    if (meta.files.size() == 1)
        arguments.bin_path.push_back(std::vector<std::string>{meta.files[0].path});
    else
    {
        for (auto & seg : meta.segments)
        {
            std::unordered_set<size_t> file_ids{};
            for (size_t seq_id : seg.seq_vec)
            {
                file_ids.emplace(meta.sequences[seq_id].file_id);
            }
            std::vector<std::string> bin_files{};
            for (size_t file_id : file_ids)
            {
                bin_files.push_back(meta.files[file_id].path);
            }
            arguments.bin_path.push_back(bin_files);
        }
    }

    // ==========================================
    // Process minimiser parameters for IBF size calculation.
    // ==========================================
    if ((parser.is_option_set("kmer-count-min") || parser.is_option_set("kmer-count-max")) && parser.is_option_set("use-filesize-dependent-cutoff"))
        throw sharg::parser_error{"You cannot use both --kmer-count-cutoff and --use-filesize-dependent-cutoff."};
        
    try
    {
        sharg::output_file_validator{sharg::output_file_open_options::open_or_create}(arguments.out_path);
        arguments.out_dir = arguments.out_path.parent_path();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    if (!parser.is_option_set("window"))
    {
        if (arguments.fast)
        {
            arguments.window_size = arguments.kmer_size + 2;
            arguments.input_is_minimiser = true;
        }
        else
            arguments.window_size = arguments.kmer_size;
    }

    if (arguments.fast)
        raptor::compute_minimiser(arguments);   // requires bin_path and out_dir
    
    // ==========================================
    // Find IBF size.
    // ==========================================
    if (parser.is_option_set("size"))
    {
        arguments.size.erase(std::remove(arguments.size.begin(), arguments.size.end(), ' '), arguments.size.end());

        size_t multiplier{};

        switch (std::tolower(arguments.size.back()))
        {
            case 't':
                multiplier = 8ull * 1024ull * 1024ull * 1024ull * 1024ull;
                break;
            case 'g':
                multiplier = 8ull * 1024ull * 1024ull * 1024ull;
                break;
            case 'm':
                multiplier = 8ull * 1024ull * 1024ull;
                break;
            case 'k':
                multiplier = 8ull * 1024ull;
                break;
            default:
                throw sharg::parser_error{"Use {k, m, g, t} to pass size. E.g., --size 8g."};
        }

        size_t size{};
        std::from_chars(arguments.size.data(), arguments.size.data() + arguments.size.size() - 1, size);
        size *= multiplier;
        arguments.bits = size / (((arguments.bins + 63) >> 6) << 6);
    }
    else 
    {
        arguments.fpr = meta.ibf_fpr;
        arguments.bits = raptor::compute_bin_size(arguments);
    }
    
    // ==========================================
    // Dispatch
    // ==========================================
    valik_build(arguments);
}

} // namespace valik::app
