#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <valik/argument_parsing/search.hpp>
#include <valik/index.hpp>
#include <valik/search/search.hpp>

namespace valik::app
{

void init_search_parser(sharg::parser & parser, search_arguments & arguments)
{
    init_shared_meta(parser);
    parser.add_option(arguments.index_file,
                      sharg::config{.short_id = '\0',
                      .long_id = "index",
                      .description = "Provide a valid path to an IBF.",
                      .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.query_file,
                      sharg::config{.short_id = '\0',
                      .long_id = "query",
                      .description = "Provide a path to the query file.",
                      .required = true,
                      .validator = sharg::input_file_validator{{"fasta", "fa", "fq", "fastq"}}});
    parser.add_flag(arguments.split_query,
                    sharg::config{.short_id = '\0',
                    .long_id = "split-query",
                    .description = "Split a long query into segments."});
    parser.add_option(arguments.out_file,
                      sharg::config{.short_id = '\0',
                      .long_id = "output",
                      .description = "Please provide a valid path to the output.",
                      .required = true,
                      .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"gff"}}});
    parser.add_option(arguments.error_rate,
                      sharg::config{.short_id = 'e',
                      .long_id = "error-rate",
                      .description = "Choose the maximum allowed error rate of a local match.",
                      .validator = sharg::arithmetic_range_validator{0.0f, 0.2f}});
    parser.add_option(arguments.pattern_size,
                      sharg::config{.short_id = '\0',
                      .long_id = "pattern",
                      .description = "Choose the size of the approximate seed."});
    parser.add_flag(arguments.fast,
                      sharg::config{.short_id = '\0',
                      .long_id = "fast",
                      .description = "Execute the search in fast mode when few false negatives can be tolerated."});
    parser.add_flag(arguments.write_time,
                    sharg::config{.short_id = '\0',
                    .long_id = "time",
                    .description = "Write runtime log file.",
                    .advanced = true});
    parser.add_flag(arguments.verbose,
                    sharg::config{.short_id = '\0',
                    .long_id = "verbose",
                    .description = "Print verbose output.",
                    .advanced = true});
    parser.add_flag(arguments.very_verbose,
                    sharg::config{.short_id = '\0',
                    .long_id = "very-verbose",
                    .description = "Print verbose output.",
                    .advanced = true});
    parser.add_option(arguments.ref_meta_path,
                    sharg::config{.short_id = '\0',
                    .long_id = "ref-meta",
                    .description = "Path to reference metadata file created by split.",
                    .validator = sharg::input_file_validator{}});
    parser.add_flag(arguments.distribute,
                    sharg::config{.short_id = '\0',
                    .long_id = "distribute",
                    .description = "Launch parallel Stellar instances each of them with independent memory."});
    parser.add_option(arguments.threads,
                    sharg::config{.short_id = '\0',
                    .long_id = "threads",
                    .description = "Choose the number of threads.",
                    .validator = positive_integer_validator{}});
    
    /////////////////////////////////////////
    // Advanced options
    /////////////////////////////////////////
    parser.add_flag(arguments.stellar_only,
                      sharg::config{.short_id = '\0',
                      .long_id = "stellar-only",
                      .description = "Do not prefilter before Stellar search.",
                      .advanced = true});
    parser.add_flag(arguments.manual_parameters,
                      sharg::config{.short_id = '\0',
                      .long_id = "without-parameter-tuning",
                      .description = "Preprocess database without setting default parameters.",
                      .advanced = true});
    parser.add_flag(arguments.keep_best_repeats,
                    sharg::config{.short_id = '\0',
                    .long_id = "keep-best-repeats",
                    .description = "Consider only high entropy regions for queries with abundant matches.",
                    .advanced = true});
    parser.add_option(arguments.best_bin_entropy_cutoff,
                      sharg::config{.short_id = '\0',
                      .long_id = "bin-entropy-cutoff",
                      .description = "For queries with abundant matches, search only highly varied reference regions. "
                                     "Increase this value to search more of the reference. "
                                     "Use with --keep-best-repeats.",
                      .advanced = true,
                      .validator = sharg::arithmetic_range_validator{0.0, 1.0}});
    parser.add_flag(arguments.keep_all_repeats,
                    sharg::config{.short_id = '\0',
                    .long_id = "keep-all-repeats",
                    .description = "Do not filter out query matches from repeat regions. This may significantly increase the runtime.",
                    .advanced = true});
    parser.add_flag(arguments.static_threshold,
                    sharg::config{.short_id = '\0',
                    .long_id = "static-threshold",
                    .description = "Do not correct threshold to avoid many spuriously matching bins.",
                    .advanced = true});
    parser.add_option(arguments.seg_count_in,
                      sharg::config{.short_id = 'n',
                      .long_id = "seg-count",
                      .description = "The suggested number of database segments that might be adjusted by the split algorithm.",
                      .validator = positive_integer_validator{false}});
    parser.add_option(arguments.threshold,
                      sharg::config{.short_id = '\0',
                      .long_id = "threshold",
                      .description = "If set, this threshold is used instead of the kmer lemma or probabilistic minimiser threshold.",
                      .advanced = true,
                      .validator = sharg::arithmetic_range_validator{1, 500}});
    parser.add_flag(arguments.cache_thresholds,
                    sharg::config{.short_id = '\0',
                    .long_id = "cache-thresholds",
                    .description = "Stores the computed thresholds with an unique name next to the index. In the next search call "
                                   "using this option, the stored thresholds are re-used.\n"
                                   "Two files are stored:\n"
                                   "\\fBthreshold_*.bin\\fP: Depends on pattern, window, kmer/shape, errors, and tau.\n"
                                   "\\fBcorrection_*.bin\\fP: Depends on pattern, window, kmer/shape, p_max, and fpr.", 
                    .advanced = true});
    parser.add_option(arguments.query_every,
                      sharg::config{.short_id = '\0',
                      .long_id = "query-every",
                      .description = "Choose how much sequential patterns overlap. "
                                     "This determines how many potential matches are skipped in prefiltering."
                                     "--query-every 1 considers all potential matches.", 
                      .advanced = true});
    parser.add_option(arguments.cart_max_capacity,
                    sharg::config{.short_id = '\0',
                    .long_id = "cart-max-capacity",
                    .description = "Number of elements to be stored in a single cart before it is sent for processing.", 
                    .advanced = true});
    parser.add_option(arguments.max_queued_carts,
                    sharg::config{.short_id = '\0',
                    .long_id = "max-queued-carts",
                    .description = "Maximal number of carts that are full and are waiting to be processed.",
                    .advanced = true});

    parser.add_subsection("Stellar options");
    parser.add_option(arguments.minLength,
                      sharg::config{.short_id = '\0',
                      .long_id = "minLength",
                      .description = "Choose the minimium length of a local alignment."});
    parser.add_option(arguments.disableThresh,
                    sharg::config{.short_id = '\0',
                    .long_id = "disableThresh",
                    .description = "STELLAR: Maximal number of verified SWIFT filter matches before disabling verification for one query sequence.",
                    .advanced = true,
                    .validator = sharg::arithmetic_range_validator{1, 10000}});
    parser.add_option(arguments.compactThresh,
                    sharg::config{.short_id = 's',
                    .long_id = "sortThresh",
                    .description = "STELLAR: Number of matches triggering removal of duplicates. Choose a smaller value for saving space.", 
                    .advanced = true});
    parser.add_option(arguments.qGram,
                    sharg::config{.short_id = 'q',
                    .long_id = "stellar-kmer",
                    .description = "STELLAR: Length of the q-grams in the SWIFT algorithm.",
                    .advanced = true,
                    .validator = sharg::arithmetic_range_validator{1, 31}});
    parser.add_option(arguments.qgramAbundanceCut,
                    sharg::config{.short_id = 'c',
                    .long_id = "abundanceCut",
                    .description = "STELLAR: k-mer overabundance cut ratio.",
                    .advanced = true,
                    .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(arguments.maxRepeatPeriod,
                    sharg::config{.short_id = '\0',
                    .long_id = "repeatPeriod",
                    .description = "STELLAR: Maximal period of low complexity repeats to be filtered.",
                    .advanced = true,
                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(arguments.minRepeatLength,
                    sharg::config{.short_id = '\0',
                    .long_id = "repeatLength",
                    .description = "STELLAR: Minimal length of low complexity repeats to be filtered.",
                    .advanced = true,
                    .validator = sharg::arithmetic_range_validator{1u, std::numeric_limits<uint32_t>::max()}});
    parser.add_option(arguments.xDrop,
                    sharg::config{.short_id = 'x',
                    .long_id = "xDrop",
                    .description = "STELLAR: Maximal x-drop for extension.", 
                    .advanced = true});
    parser.add_option(arguments.strVerificationMethod,
                    sharg::config{.short_id = '\0',
                    .long_id = "verification",
                    .description = "STELLAR: Verification strategy: exact or bestLocal.",
                    .advanced = true,
                    .validator = sharg::value_list_validator{"exact", "bestLocal"}});
    parser.add_option(arguments.numMatches,
                    sharg::config{.short_id = '\0',
                    .long_id = "numMatches",
                    .description = "STELLAR: Maximal number of kept matches per query and database." 
                                   "If STELLAR finds more matches, only the longest ones are kept.",
                    .advanced = true});
}

void run_search(sharg::parser & parser)
{
    search_arguments arguments{};

    init_search_parser(parser, arguments);

    try_parsing(parser);

    if (arguments.very_verbose)
        arguments.verbose = true;

    if (arguments.manual_parameters && !parser.is_option_set("ref-meta") && arguments.split_query)
    {
        if (!parser.is_option_set("pattern") || !parser.is_option_set("seg-count"))
            throw std::runtime_error("Provide --ref-meta to deduce parameter values or provide --seg-count and --pattern manually.");
    }

    // ==========================================
    // Process --seg-count.
    // ==========================================
    if (parser.is_option_set("seg-count"))
    {
        arguments.split_query = true;
    }

    // ==========================================
    // Checking advanced Stellar parameters.
    // ==========================================
    sharg::input_file_validator{}(arguments.query_file);
    if (parser.is_option_set("disableThresh") && parser.is_option_set("numMatches"))
    {
        if (arguments.numMatches > arguments.disableThresh)
            throw sharg::validation_error{"Disabling verification for queries with disableThresh=" + std::to_string(arguments.disableThresh) + 
                                          " or more matches.\n The threshold for maximum numer of matches per query numMatches=" + std::to_string(arguments.numMatches) + 
                                          " can not be larger."};
    }

    if (arguments.numMatches > arguments.compactThresh)
        throw sharg::validation_error{"Invalid parameter values: Please choose numMatches <= sortThresh.\n"};


    if (!arguments.stellar_only && !parser.is_option_set("index"))
    {
        throw sharg::parser_error{"Option --index is required but not set."};
    }

    if (!arguments.stellar_only)
    {
        // ==========================================
        // Read window and kmer size, and the bin paths.
        // ==========================================
        std::ifstream is{arguments.index_file.string(),std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        valik_index<> tmp{};
        tmp.load_parameters(iarchive);
        arguments.shape = tmp.shape();
        arguments.shape_size = arguments.shape.size();
        arguments.shape_weight = arguments.shape.count();
        arguments.window_size = tmp.window_size();
        arguments.bin_path = tmp.bin_path();
        if (arguments.bin_path.size() > 1)
            arguments.distribute = true;
    }

    // ==========================================
    // Process --pattern.
    // ==========================================
    if (parser.is_option_set("pattern"))
    {
        if (!arguments.manual_parameters)
        {
            std::cerr << "WARNING: pattern size (minimum match length) will be adjusted to match database metadata. "
                      << "Set --without-parameter-tuning to force manual input.\n"; 
        }
        else
        {
            if (arguments.pattern_size < arguments.window_size)
                throw sharg::validation_error{"The minimiser window cannot be bigger than the pattern."};
        }
    }
    else
    {
        if (arguments.manual_parameters)
            throw sharg::validation_error{"Input --pattern size or deduce parameter by providing --ref-meta."};
    }

    arguments.errors = std::ceil(arguments.error_rate * arguments.pattern_size);
    
    // ==========================================
    // Process manual threshold.
    // ==========================================
    if (parser.is_option_set("threshold") && !arguments.stellar_only)
    {
        arguments.threshold_was_set = true;
        arguments.threshold_percentage = arguments.threshold / (double) (arguments.pattern_size - arguments.shape.size() + 1);
    }

    // ==========================================
    // Process reference metadata path.
    // ==========================================
    if (parser.is_option_set("ref-meta"))
    {
        // Create temporary file path for merging parallel Stellar runs.
        arguments.all_matches = arguments.out_file;
        arguments.all_matches += ".preliminary";
    }
    else
    {
        if (arguments.split_query && !arguments.manual_parameters)
        {
            throw std::runtime_error("Provide --ref-meta to deduce parameter values or provide --seg-count and --pattern manually.");
        }
    }

    if (!arguments.manual_parameters)
    {
        // ==========================================
        // Extract data from reference metadata.
        // ==========================================
        if (!parser.is_option_set("ref-meta"))
            throw sharg::validation_error("Provide --ref-meta to deduce suitable search parameters or set --without-parameter-tuning and --pattern size.");

        std::filesystem::path search_profile_file{arguments.ref_meta_path};
        search_profile_file.replace_extension("arg");
        sharg::input_file_validator argument_input_validator{{"arg"}};
        argument_input_validator(search_profile_file);
        search_kmer_profile search_profile{search_profile_file};

        arguments.pattern_size = search_profile.l;
        arguments.errors = std::ceil(arguments.error_rate * arguments.pattern_size);    // update based on pattern size in metadata
        search_error_profile error_profile = search_profile.get_error_profile(arguments.errors);
        // seg_count is inferred in metadata constructor
        arguments.search_type = error_profile.search_type;
        if (parser.is_option_set("threshold"))
        {
            utilities::kmer kmer{arguments.shape};
            auto lemma_thresh = kmer.lemma_threshold(arguments.pattern_size, arguments.errors);
            if (arguments.threshold > lemma_thresh)
                arguments.search_type = search_kind::HEURISTIC;
            else 
            {
                arguments.search_type = search_kind::LEMMA;
                if (arguments.threshold < lemma_thresh)
                    std::cerr << "[Warning] The chosen threshold is less than the k-mer lemma threshold. Ignore this warning if this was deliberate.";
            }
        }
        if (arguments.stellar_only)
            arguments.search_type = search_kind::STELLAR;

        if (arguments.search_type == search_kind::STELLAR)
        {
            std::cout << "Can not prefilter matches of length " << std::to_string(error_profile.pattern.l) << 
                                     " with " << std::to_string(error_profile.pattern.e) << " errors. Searching without prefiltering.\n";
            arguments.fnr = 0.0;
        }
        else
        {
            arguments.max_segment_len = max_segment_len(error_profile.fp_per_pattern, arguments.pattern_size, arguments.query_every);
            if (!parser.is_option_set("threshold"))
                arguments.threshold = error_profile.params.t;
            arguments.threshold_was_set = true;  // use raptor::threshold_kinds::percentage
            if (arguments.threshold > arguments.pattern_size - arguments.shape.size() + 1)
                throw sharg::validation_error("Threshold can not be larger than the number of k-mers per pattern.");
            arguments.threshold_percentage = arguments.threshold / (double) (arguments.pattern_size - arguments.shape.size() + 1);
            arguments.fnr = error_profile.fnr;

            if (arguments.window_size > arguments.shape_size)
                arguments.search_type = search_kind::MINIMISER;
        }

        // ==========================================
        // More checks.
        // ==========================================
        if (parser.is_option_set("query-every"))
        {
            if (arguments.query_every > arguments.pattern_size)
                throw sharg::validation_error("Reduce --query-every so that all positions in the query would be considered at least once."); 
        }
        else
        {
            if (arguments.fast)
                arguments.query_every = 3;
        }

        // ==========================================
        // Set strict thresholding parameters for fast mode.
        // ==========================================
        if (arguments.fast)
            arguments.tau = 0.99999;

        if (arguments.fast)
            arguments.p_max = 0.05;
    }

    if (!parser.is_option_set("minLength"))
        arguments.minLength = arguments.pattern_size;

    // ==========================================
    // Dispatch
    // ==========================================
    valik_search(arguments);
}

} // namespace valik::app
