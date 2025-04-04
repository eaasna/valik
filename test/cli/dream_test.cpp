#include <filesystem>
#include <fstream>
#include <ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include "cli_test.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// DREAM short shared memory ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(dream_short_search, short_shared_mem)
{
    auto const [number_of_errors] = GetParam();
    size_t pattern_size = 50;
    float error_rate = (float) number_of_errors / (float) pattern_size;
    float max_error_rate = 0.04;

    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);

    std::filesystem::path ref_meta_path = "ref.bin";
    std::filesystem::path index_path = "ref.ibf";

    cli_test_result const build = execute_app("valik", "build",
                                                       data("ref.fasta"),
                                                       "--output ", index_path, 
                                                       "--fpr 0.001", 
                                                       "--pattern ", std::to_string(pattern_size),
                                                       "--error-rate ", std::to_string(max_error_rate));
    EXPECT_EQ(build.exit_code, 0);
    valik::metadata reference(ref_meta_path);

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", index_path,
                                                        "--query ", data("query.fasta"),
                                                        "--repeatPeriod 1",
                                                        "--repeatLength 10", 
                                                        "--numMatches 2");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 4 database sequences.\n"});
    EXPECT_EQ(result.err, std::string{});

    auto distributed = valik::read_alignment_output<valik::stellar_match>(search_result_path(number_of_errors), reference, std::ios::binary);
    auto local = valik::read_alignment_output<valik::stellar_match>("search.gff", reference);

    compare_gff_out(distributed, local);
}

INSTANTIATE_TEST_SUITE_P(short_shared_memory_dream_suite,
                         dream_short_search,
                         testing::Combine(testing::Values(1, 2)),
                         [] (testing::TestParamInfo<dream_short_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_error";
                             return name;
                         });

TEST_F(dream_short_search, no_matches)
{
    setup_tmp_dir();
    size_t pattern_size = 50;

    setenv("VALIK_MERGE", "cat", true);

    std::filesystem::path ref_meta_path = "ref.bin";
    std::filesystem::path index_path = "ref.ibf";
    
    cli_test_result const build = execute_app("valik", "build",
                                                       data("ref.fasta"),
                                                       "--pattern ", std::to_string(pattern_size),
                                                       "--output ", index_path);
    EXPECT_EQ(build.exit_code, 0);
    EXPECT_EQ(build.err, std::string{});

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--error-rate 0",
                                                        "--index ", index_path,
                                                        "--query ", data("dummy_reads.fastq"));

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 4 database sequences.\n"});
    EXPECT_EQ(result.err, std::string{}); // Stellar shortens query IDs until the first whitespace

    auto actual = string_list_from_file("search.gff");

    EXPECT_EQ(0, actual.size());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// DREAM split shared memory ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(dream_split_search, split_shared_mem)
{
    auto const [number_of_errors] = GetParam();
    size_t pattern_size = 50;
    float error_rate = (float) number_of_errors / (float) pattern_size;
    float max_error_rate = 0.04;

    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);

    std::filesystem::path ref_meta_path = "ref.bin";
    std::filesystem::path index_path = "ref.ibf";

    cli_test_result const build = execute_app("valik", "build",
                                                       data("ref.fasta"),
                                                       "--fpr 0.001",
                                                       "--pattern ", std::to_string(pattern_size), 
                                                       "--error-rate ", std::to_string(max_error_rate),
                                                       "--output ", index_path);

    EXPECT_EQ(build.exit_code, 0);
    valik::metadata reference(ref_meta_path);

    cli_test_result const search = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--split-query",
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", index_path,
                                                        "--query ", data("query.fasta"),
                                                        "--repeatPeriod 1",
                                                        "--repeatLength 10",
                                                        "--numMatches 2");

    EXPECT_EQ(search.exit_code, 0);
    EXPECT_EQ(search.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 4 database sequences.\n"});
    EXPECT_EQ(search.err, std::string{});

    auto distributed = valik::read_alignment_output<valik::stellar_match>(search_result_path(number_of_errors), reference, std::ios::binary);
    auto local = valik::read_alignment_output<valik::stellar_match>("search.gff", reference);
    
    compare_gff_out(distributed, local);
}

INSTANTIATE_TEST_SUITE_P(split_shared_memory_dream_suite,
                         dream_split_search,
                         testing::Combine(testing::Values(1, 2)),
                         [] (testing::TestParamInfo<dream_split_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_error";
                             return name;
                         });


TEST_P(dream_split_search, split_shaped_kmer)
{
    auto const [number_of_errors] = GetParam();
    size_t pattern_size = 50;
    float error_rate = (float) number_of_errors / (float) pattern_size;
    float max_error_rate = 0.04;

    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);

    std::filesystem::path ref_meta_path = "ref.bin";
    std::filesystem::path index_path = "ref.ibf";

    cli_test_result const build = execute_app("valik", "build",
                                                       data("ref.fasta"),
                                                       "--output ", index_path,
                                                       "--fpr 0.001",
                                                       "--shape 1111110110110111111",
                                                       "--pattern ", std::to_string(pattern_size), 
                                                       "--error-rate ", std::to_string(max_error_rate));

    EXPECT_EQ(build.exit_code, 0);
    valik::metadata reference(ref_meta_path);

    cli_test_result const search = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--split-query",
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", index_path,
                                                        "--query ", data("query.fasta"),
                                                        "--repeatPeriod 1",
                                                        "--repeatLength 10", 
                                                        "--numMatches 2");

    EXPECT_EQ(search.exit_code, 0);
    EXPECT_EQ(search.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 4 database sequences.\n"});
    EXPECT_EQ(search.err, std::string{});

    auto distributed = valik::read_alignment_output<valik::stellar_match>(search_result_path(number_of_errors), reference, std::ios::binary);
    auto local = valik::read_alignment_output<valik::stellar_match>("search.gff", reference);
    
    compare_gff_out(distributed, local);
}

INSTANTIATE_TEST_SUITE_P(split_shared_memory_gapped_dream_suite,
                         dream_split_search,
                         testing::Combine(testing::Values(1, 2)),
                         [] (testing::TestParamInfo<dream_split_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_error";
                             return name;
                         });
