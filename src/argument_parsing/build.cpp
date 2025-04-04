#include <valik/split/split.hpp>
#include <valik/argument_parsing/build.hpp>

namespace valik::app
{

void init_build_parser(sharg::parser & parser, build_arguments & arguments)
{
    param_space space{};
    init_shared_meta(parser);

    /////////////////////////////////////////
    // Split options
    /////////////////////////////////////////
    parser.add_positional_option(arguments.db_file,
                      sharg::config{.description = "File containing database sequences. If splitting --metagenome provide a text file with a list of cluster paths.",
                      .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.pattern_size,
                      sharg::config{.short_id = '\0',
                      .long_id = "pattern",
                      .description = "Choose how much consecutive segments overlap. This is the minimum length of a local alignment.",
                      .validator = sharg::arithmetic_range_validator{10u, 1000u}});
    parser.add_option(arguments.error_rate,
                      sharg::config{.short_id = 'e',
                      .long_id = "error-rate",
                      .description = "Choose the upper bound for the maximum allowed error rate of a local match.",
                      .validator = sharg::arithmetic_range_validator{0.0f, 0.2f}});
    parser.add_option(arguments.fpr,
                      sharg::config{.short_id = '\0',
                      .long_id = "fpr",
                      .description = "False positive rate of IBF.",
                      .validator = sharg::arithmetic_range_validator{0.0001f, 0.5f}});
    parser.add_option(arguments.kmer_size,
                      sharg::config{.short_id = 'k',
                      .long_id = "kmer",
                      .description = "Choose the ungapped kmer size.",
                      .validator = sharg::arithmetic_range_validator{space.min_k(), space.max_k()}});
    parser.add_option(arguments.shape_str,
                      sharg::config{.short_id = 's',
                      .long_id = "shape",
                      .description = "Choose the kmer shape. E.g 15-mer with two gaps: 111110111011111.",
                      .validator = sharg::regex_validator{"[01]+"}});
    parser.add_option(arguments.seg_count_in,
                      sharg::config{.short_id = 'n',
                      .long_id = "seg-count",
                      .description = "The suggested number of database segments that might be adjusted by the split algorithm.",
                      .validator = positive_integer_validator{false}});
    parser.add_flag(arguments.metagenome,
                      sharg::config{.short_id = '\0',
                      .long_id = "metagenome",
                      .description = "Split a clustered metagenome database. Reference input is a list of cluster paths"});

    parser.add_option(arguments.out_path,
                      sharg::config{.short_id = 'o',
                      .long_id = "output",
                      .description = "Provide an output filepath."});
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
    parser.add_flag(arguments.write_out,
                      sharg::config{.short_id = '\0',
                      .long_id = "write-out",
                      .description = "Write an output FASTA file for each reference segment or write all query segments into a single output FASTA file.",
                      .advanced = true});
    parser.add_option(arguments.window_size,
                      sharg::config{.short_id = 'w',
                      .long_id = "window",
                      .description = "Choose the window size.",
                      .advanced = true,
                      .validator = positive_integer_validator{}});
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
}

void run_build(sharg::parser & parser)
{
    build_arguments arguments{};
    init_build_parser(parser, arguments);
    try_parsing(parser);

    if (parser.is_option_set("kmer"))
    {
        if (parser.is_option_set("shape"))
        {
            throw sharg::parser_error{"Arguments --kmer and --shape are mutually exclusive."};
        }
        else
        {
            arguments.shape = seqan3::shape{seqan3::ungapped(arguments.kmer_size)};
            arguments.shape_str = std::string(arguments.kmer_size, '1');
            arguments.shape_weight = arguments.kmer_size;
        }
    }

    if (parser.is_option_set("shape"))
    {
        uint64_t bin_shape{};
        std::from_chars(arguments.shape_str.data(), arguments.shape_str.data() + arguments.shape_str.size(), bin_shape, 2);
        arguments.shape = seqan3::shape(seqan3::bin_literal{bin_shape});
        arguments.kmer_size = arguments.shape.size();
        arguments.shape_weight = arguments.shape.count();
    }

    arguments.errors = std::ceil(arguments.error_rate * arguments.pattern_size);
    // ==========================================
    // Process bin_path:
    // if building from clustered sequences each line in input corresponds to a bin
    // if building from overlapping segments all sequences in one reference file
    // ==========================================
    auto sequence_file_validator{bin_validator{}.sequence_file_validator};

    if (arguments.metagenome)
    {
        std::ifstream istrm{arguments.db_file};
        std::string bin_path;
        while (std::getline(istrm, bin_path))
        {
            if (!bin_path.empty())
            {
                sequence_file_validator(bin_path);
                arguments.bin_path.emplace_back(bin_path);
            }
        }
    }
    else
    {
        sequence_file_validator(arguments.db_file);
        std::string bin_string{arguments.db_file.string()};
        arguments.bin_path.emplace_back(bin_string);
    }

    if (!parser.is_option_set("output"))
    {
        arguments.out_path = arguments.db_file;
        arguments.out_path.replace_extension("index");
    }
    sharg::output_file_validator{sharg::output_file_open_options::open_or_create}(arguments.out_path);
    
    arguments.ref_meta_path = arguments.out_path;
    arguments.ref_meta_path.replace_extension("bin");

    if (parser.is_option_set("seg-count") && !arguments.only_split)
    {
        std::cerr << "WARNING: seg count will be adjusted to the next multiple of 64. "
                  << "Set --without-parameter-tuning to force manual input.\n";
    }

    // ==========================================
    // Dispatch split
    // ==========================================
    auto meta = valik_split(arguments);
    arguments.bins = meta.seg_count;

    if (!arguments.manual_parameters)
    {
        std::filesystem::path search_profile_file{arguments.ref_meta_path};
        search_profile_file.replace_extension("arg");
        sharg::input_file_validator argument_input_validator{{"arg"}};
        argument_input_validator(search_profile_file);
        search_kmer_profile search_profile{search_profile_file};

        if (parser.is_option_set("kmer"))
        {
            seqan3::debug_stream << "WARNING: kmer size k=" << arguments.kmer_size << " will be updated to " << search_profile.kmer.size() 
                                 << ". Set --without-parameter-tuning to force manual input.";
        }
        if (parser.is_option_set("window"))
        {
            seqan3::debug_stream << "WARNING: window size w=" << arguments.window_size << " will be updated. " 
                                 << "Set --without-parameter-tuning to force manual input.";
        }
        
        arguments.kmer_size = search_profile.kmer.size();
        arguments.window_size = search_profile.kmer.size();
        arguments.shape = search_profile.kmer.shape;
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

    try
    {
        arguments.out_dir = arguments.out_path.parent_path();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    // ==========================================
    // Process minimiser parameters for IBF size calculation.
    // ==========================================
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
