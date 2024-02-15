#include <seqan3/core/debug_stream.hpp>

#include <valik/argument_parsing/split.hpp>
#include <valik/split/split.hpp>

namespace valik::app
{

void init_split_parser(sharg::parser & parser, split_arguments & arguments)
{
    init_shared_meta(parser);
    parser.add_positional_option(arguments.db_file,
                      sharg::config{.description = "File containing database sequences. If splitting --metagenome provide a text file with a list of cluster paths.",
                      .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.meta_out,
                      sharg::config{.short_id = 'o',
                      .long_id = "out",
                      .description = "Please provide a valid path to the database metadata output.",
                      .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create}});
    parser.add_option(arguments.pattern_size,
                      sharg::config{.short_id = '\0',
                      .long_id = "pattern",
                      .description = "Choose how much consecutive segments overlap. This is the minimum length of a local alignment.",
                      .validator = positive_integer_validator{true}});
    parser.add_option(arguments.error_rate,
                      sharg::config{.short_id = 'e',
                      .long_id = "error-rate",
                      .description = "Choose the upper bound for the maximum allowed error rate of a local match.",
                      .validator = sharg::arithmetic_range_validator{0.0f, 0.1f}});
    parser.add_option(arguments.kmer_size,
                      sharg::config{.short_id = 'k',
                      .long_id = "kmer",
                      .description = "Choose the kmer size.",
                      .validator = positive_integer_validator{true}});
    parser.add_option(arguments.seg_count_in,
                      sharg::config{.short_id = 'n',
                      .long_id = "seg-count",
                      .description = "The suggested number of database segments that might be adjusted by the split algorithm.",
                      .validator = positive_integer_validator{false}});
    parser.add_flag(arguments.split_index,
                      sharg::config{.short_id = '\0',
                      .long_id = "split-index",
                      .description = "Adjust the suggested segment count to create a multiple of 64 segments instead. This is suitable for building an IBF."});
    parser.add_flag(arguments.metagenome,
                      sharg::config{.short_id = '\0',
                      .long_id = "metagenome",
                      .description = "Split a clustered metagenome database. Reference input is a list of cluster paths"});
    parser.add_option(arguments.ref_meta_path,
                    sharg::config{.short_id = '\0',
                    .long_id = "ref-meta",
                    .description = "Path to reference metadata file created by split.",
                    .validator = sharg::input_file_validator{}});
    parser.add_flag(arguments.write_out,
                      sharg::config{.short_id = '\0',
                      .long_id = "write-out",
                      .description = "Write an output FASTA file for each reference segment or write all query segments into a single output FASTA file..",
                      .advanced = true});
    parser.add_flag(arguments.only_split,
                      sharg::config{.short_id = '\0',
                      .long_id = "without-parameter-tuning",
                      .description = "Preprocess database without setting default parameters.",
                      .advanced = true});                      
    parser.add_flag(arguments.verbose,
                      sharg::config{.short_id = 'v',
                      .long_id = "verbose",
                      .description = "Print verbose output.",
                      .advanced = true});

}

void run_split(sharg::parser & parser)
{
    split_arguments arguments{};
    init_split_parser(parser, arguments);
    try_parsing(parser);

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
        std::string line;
        while (std::getline(istrm, line))
        {
            if (!line.empty())
            {
                std::stringstream linestrm(line);
                std::string filepath;
                std::vector<std::string> bin_paths;
                while (std::getline(linestrm, filepath, '\t'))
                {
                    sequence_file_validator(filepath);
                    bin_paths.emplace_back(filepath);
                }

                arguments.bin_path.emplace_back(bin_paths);
            }
        }
    }
    else
    {
        sequence_file_validator(arguments.db_file);
        std::string bin_string{arguments.db_file.string()};
        arguments.bin_path.emplace_back(std::vector<std::string>{bin_string});
    }

    if (!parser.is_option_set("out"))
    {
        arguments.meta_out = arguments.db_file;
        arguments.meta_out.replace_extension("bin");
    }

    if (arguments.metagenome)
        arguments.split_index = true;

    if (!arguments.split_index && !arguments.only_split && !parser.is_option_set("ref-meta"))
        throw sharg::parser_error{"Need to provide path to reference metadata to process a query database."};

    // ==========================================
    // Dispatch
    // ==========================================
    valik_split(arguments);
}

} // namespace valik::app
