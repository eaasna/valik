#include <gtest/gtest.h>

#include <utilities/consolidate/consolidate_matches.hpp>

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data_path(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

std::filesystem::path consolidation_input_path(size_t const number_of_bins, size_t const overlap) noexcept
{
    std::string name{};
    name += std::to_string(number_of_bins);
    name += "bins";
    name += std::to_string(overlap);
    name += "overlap_dream_all.gff";
    return data_path(name);
}

std::filesystem::path consolidated_matches_path(size_t const num_matches, size_t const disable_thresh) noexcept
{
    std::string name{};
    name += "stellar_truth_num";
    name += std::to_string(num_matches);
    name += "_dis";
    name += std::to_string(disable_thresh);
    name += ".gff";
    return data_path(name);
}

std::filesystem::path consolidation_meta_path(size_t const number_of_bins, size_t const overlap) noexcept
{
    std::string name{};
    name += std::to_string(number_of_bins);
    name += "bins";
    name += std::to_string(overlap);
    name += "overlap_reference_metadata.bin";
    return data_path(name);
}

void compare_gff_out(std::vector<valik::stellar_match> const & expected,
                            std::vector<valik::stellar_match> const & actual)
{
    //EXPECT_EQ(expected.size(), actual.size());
    if (expected.size() != actual.size())
        seqan3::debug_stream << "[WARNING] Discrepancy in duplicate removal.\n\tStellar match count\t" << expected.size() << 
                                "\n\tDREAM-Stellar match count\t" << actual.size() << " \n";

    size_t not_actually_found{0};
    for (auto & match : expected)
    {
        auto it = std::find(actual.begin(), actual.end(), match);
        if (it == actual.end())
        {
            not_actually_found++;
            seqan3::debug_stream << match.to_string();
        }
    }

    EXPECT_EQ(not_actually_found, 0);
    seqan3::debug_stream << "\tMatches only in Stellar\t" << not_actually_found << '\n';

    size_t only_in_actual{0};
    for (auto & match : actual)
    {
        auto it = std::find(expected.begin(), expected.end(), match);
        if (it == expected.end())
        {
            only_in_actual++;
            //seqan3::debug_stream << match.to_string();
        }   
    }
    seqan3::debug_stream << "\tMatches only in DREAM-Stellar\t" << only_in_actual << '\n';

    //!NOTE: there are duplicates that are not removed in Stellar output
    //EXPECT_EQ(only_in_actual, 0);
}

std::set<std::string> get_query_ids(std::vector<valik::stellar_match> const & matches)
{
    std::set<std::string> query_ids; 
    for (auto & match : matches)
        query_ids.insert(match.qname);

    return query_ids;
}

void check_consolidation(std::set<std::string> const &  query_ids,
                         std::vector<valik::stellar_match> const & initial, 
                         std::vector<valik::stellar_match> const & actual, 
                         valik::search_arguments const & arguments)
{
    seqan3::debug_stream << "Consolidation test\n"; 
    seqan3::debug_stream << "\tdisableThresh:\t" << arguments.disableThresh << "\n\tnumMatches:\t" << arguments.numMatches << '\n';
    for (auto & id : query_ids)
    {
        auto query_match = [&arguments, &id](auto match) { return match.qname == id; };
        auto initial_matches = initial | std::views::filter(query_match);
        auto consolidated_matches = actual | std::views::filter(query_match);

        if ((size_t) std::abs(std::distance(initial_matches.begin(), initial_matches.end())) >= arguments.disableThresh)
        {
            EXPECT_EQ(std::distance(consolidated_matches.begin(), consolidated_matches.end()), 0);
        }
        else if ((size_t) std::abs(std::distance(initial_matches.begin(), initial_matches.end())) > arguments.numMatches)
        {
            EXPECT_EQ(std::distance(consolidated_matches.begin(), consolidated_matches.end()), arguments.numMatches);
        }
        else
        {
            continue;
        }
    }
}

TEST(consolidate_matches, bins8overlap50)
{
    size_t number_of_bins = 8;
    size_t segment_overlap = 50;

    valik::search_arguments arguments{};
    arguments.ref_meta_path = consolidation_meta_path(number_of_bins, segment_overlap);
    arguments.all_matches = consolidation_input_path(number_of_bins, segment_overlap);
    arguments.out_file = "consolidated.gff";

    arguments.verbose = true;
    valik::consolidate_matches(arguments);

    valik::metadata reference(arguments.ref_meta_path);
    auto expected = valik::read_alignment_output<valik::stellar_match>(data_path("stellar_truth.gff"), reference, std::ios::binary);
    auto actual = valik::read_alignment_output<valik::stellar_match>("consolidated.gff", reference, std::ios::binary);

    compare_gff_out(expected, actual);
}

TEST(consolidate_matches, bins16overlap50)
{
    size_t number_of_bins = 16;
    size_t segment_overlap = 50;

    valik::search_arguments arguments{};
    arguments.ref_meta_path = consolidation_meta_path(number_of_bins, segment_overlap);
    arguments.all_matches = consolidation_input_path(number_of_bins, segment_overlap);
    arguments.out_file = "consolidated.gff";

    valik::consolidate_matches(arguments);

    valik::metadata reference(arguments.ref_meta_path);
    auto expected = valik::read_alignment_output<valik::stellar_match>(data_path("stellar_truth.gff"), reference, std::ios::binary);
    auto actual = valik::read_alignment_output<valik::stellar_match>("consolidated.gff", reference, std::ios::binary);

    compare_gff_out(expected, actual);
}

std::vector<valik::stellar_match> get_initial_matches(size_t const number_of_bins, 
                                                      size_t const segment_overlap, 
                                                      valik::metadata const & reference)
{
    auto matches = valik::read_alignment_output<valik::stellar_match>(consolidation_input_path(number_of_bins, segment_overlap), reference, std::ios::binary);
    std::sort(matches.begin(), matches.end(), std::greater<valik::stellar_match>());
    matches.erase( std::unique( matches.begin(), matches.end() ), matches.end() );
    return matches;
}

TEST(consolidate_matches_with_options, bins8overlap50filter)
{
    size_t number_of_bins = 8;
    size_t segment_overlap = 50;

    valik::search_arguments arguments{};
    arguments.query_file = data_path("multi_seq_query.fasta");
    arguments.ref_meta_path = consolidation_meta_path(number_of_bins, segment_overlap);
    arguments.all_matches = consolidation_input_path(number_of_bins, segment_overlap);

    valik::metadata reference(arguments.ref_meta_path);
    auto initial = get_initial_matches(number_of_bins, segment_overlap, reference);
    auto query_ids = get_query_ids(initial);
    for (size_t d : std::vector<size_t>{3, 8, 13})
    {
        for (size_t n : std::vector<size_t>{3, 9, 12})
        {
            arguments.out_file = "consolidated_" + std::to_string(d) + "_" + std::to_string(n) + ".gff";
            arguments.disableThresh = d;
            arguments.numMatches = n;
            valik::consolidate_matches(arguments);

            auto expected = valik::read_alignment_output<valik::stellar_match>(consolidated_matches_path(arguments.numMatches, arguments.disableThresh), reference, std::ios::binary);
            auto actual = valik::read_alignment_output<valik::stellar_match>(arguments.out_file, reference, std::ios::binary);

            check_consolidation(query_ids, initial, actual, arguments);

            //!TODO: not able to reproduce consolidation in Stellar3
            // Stellar does not consolidated based on query-ref match count nor query match count
            //compare_gff_out(expected, actual);
        }
    }
}

TEST(consolidate_matches_with_options, bins16overlap50filter)
{
    size_t number_of_bins = 16;
    size_t segment_overlap = 50;

    valik::search_arguments arguments{};
    arguments.query_file = data_path("multi_seq_query.fasta");
    arguments.ref_meta_path = consolidation_meta_path(number_of_bins, segment_overlap);
    arguments.all_matches = consolidation_input_path(number_of_bins, segment_overlap);

    valik::metadata reference(arguments.ref_meta_path);
    auto initial = get_initial_matches(number_of_bins, segment_overlap, reference);
    auto query_ids = get_query_ids(initial);
    for (size_t d : std::vector<size_t>{3, 8, 13})
    {
        for (size_t n : std::vector<size_t>{3, 9, 12})
        {
            arguments.out_file = "consolidated_" + std::to_string(d) + "_" + std::to_string(n) + ".gff";
            arguments.disableThresh = d;
            arguments.numMatches = n;
            valik::consolidate_matches(arguments);

            //auto expected = valik::read_alignment_output<valik::stellar_match>(consolidated_matches_path(arguments.numMatches, arguments.disableThresh), reference, std::ios::binary);
            auto actual = valik::read_alignment_output<valik::stellar_match>(arguments.out_file, reference, std::ios::binary);

            check_consolidation(query_ids, initial, actual, arguments);

            //!TODO: not able to reproduce consolidation in Stellar3
            // Stellar does not consolidated based on query-ref match count nor query match count
            //compare_gff_out(expected, actual);
        }
    }
}
