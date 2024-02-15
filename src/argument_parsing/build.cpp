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
    parser.add_option(arguments.size,
                      sharg::config{.short_id = '\0',
                      .long_id = "size",
                      .description = "Choose the size of the resulting IBF.",
                      .required = true,
                      .validator = size_validator{"\\d+\\s{0,1}[k,m,g,t,K,M,G,T]"}});
    parser.add_flag(arguments.compressed,
                    sharg::config{.short_id = '\0',
                    .long_id = "compressed",
                    .description = "Build a compressed IBF."});
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
    
    /////////////////////////////////////////
    // Advanced options
    /////////////////////////////////////////
    parser.add_option(arguments.window_size,
                      sharg::config{.short_id = '\0',
                      .long_id = "window",
                      .description = "Choose the window size.",
                      .advanced = true,
                      .validator = positive_integer_validator{}});
    parser.add_option(arguments.kmer_size,
                      sharg::config{.short_id = '\0',
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
}

void run_build(sharg::parser & parser)
{
    build_arguments arguments{};
    init_build_parser(parser, arguments);
    try_parsing(parser);

    if (!parser.is_option_set("kmer") && !parser.is_option_set("window"))
    {
        //!TODO: read in parameter metadata file
        arguments.kmer_size = 12;
        arguments.window_size = 12;
    }

    // ==========================================
    // Various checks.
    // ==========================================
    arguments.shape = seqan3::shape{seqan3::ungapped{arguments.kmer_size}};
    arguments.shape_weight = arguments.shape.count();

    if (parser.is_option_set("window"))
    {
        if (arguments.kmer_size > arguments.window_size)
            throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};
    }
    else
    {
        if (arguments.fast)
            arguments.window_size = arguments.shape.size() + 2;
        else
            arguments.window_size = arguments.shape.size();
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

    try
    {
        sharg::output_file_validator{sharg::output_file_open_options::open_or_create}(arguments.out_path);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    // ==========================================
    // Process --size.
    // ==========================================
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

    // ==========================================
    // Dispatch
    // ==========================================
    valik_build(arguments);
}

} // namespace valik::app
