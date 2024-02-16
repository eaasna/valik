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
    auto const [number_of_bins, window_size, number_of_errors] = GetParam();
    size_t pattern_size = 50;
    float error_rate = (float) number_of_errors / (float) pattern_size;

    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);

    std::filesystem::path ref_meta_path = "ref_meta.bin";
    std::filesystem::path index_path = ibf_path(number_of_bins, window_size);
    
    cli_test_result const split_ref = execute_app("valik", "split",
                                                        data("ref.fasta"),
                                                        "--out", ref_meta_path,
                                                        "--seg-count ", std::to_string(number_of_bins),
                                                        "--pattern 150", 
                                                        "--without-parameter-tuning");
    EXPECT_EQ(split_ref.exit_code, 0);
    valik::metadata reference(ref_meta_path);

    cli_test_result const build = execute_app("valik", "build",
                                              "--window", std::to_string(window_size),
                                              "--kmer 13",
                                              "--size 32k",
                                              "--ref-meta", ref_meta_path,
                                              "--output ", index_path);
    EXPECT_EQ(build.exit_code, 0);

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--pattern ", std::to_string(pattern_size),
                                                        "--overlap ", std::to_string(pattern_size - 1),
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", index_path,
                                                        "--query ", data("query.fasta"),
                                                        "--ref-meta", ref_meta_path, 
                                                        "--repeatPeriod 1",
                                                        "--repeatLength 10", 
                                                        "--numMatches 2");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 4 database sequences.\n"});
    EXPECT_EQ(result.err, std::string{});

    auto distributed = valik::read_stellar_output(search_result_path(number_of_bins, window_size, number_of_errors), reference, std::ios::binary);
    auto local = valik::read_stellar_output("search.gff", reference);

    compare_gff_out(distributed, local);
}

INSTANTIATE_TEST_SUITE_P(short_shared_memory_dream_suite,
                         dream_short_search,
                         testing::Combine(testing::Values(4, 16), testing::Values(13, 15), testing::Values(1)),
                         [] (testing::TestParamInfo<dream_short_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });

TEST_F(dream_short_search, no_matches)
{
    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);

    size_t number_of_bins = 4;
    size_t window_size = 15;
    std::filesystem::path ref_meta_path = "ref_meta.bin";
    std::filesystem::path index_path = ibf_path(number_of_bins, window_size);
    
    cli_test_result const split_ref = execute_app("valik", "split",
                                                        data("ref.fasta"),
                                                        "--out", ref_meta_path,
                                                        "--seg-count ", std::to_string(number_of_bins),
                                                        "--pattern 150", 
                                                        "--without-parameter-tuning");
    EXPECT_EQ(split_ref.exit_code, 0);

    cli_test_result const build = execute_app("valik", "build",
                                                       "--window ", std::to_string(window_size),
                                                       "--kmer 13",
                                                       "--size 32k",
                                                       "--ref-meta", ref_meta_path,
                                                       "--output ", index_path);
    EXPECT_EQ(build.exit_code, 0);
    EXPECT_EQ(build.err, std::string{});

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--pattern 50",
                                                        "--overlap 49",
                                                        "--error-rate 0",
                                                        "--index ", ibf_path(number_of_bins, window_size),
                                                        "--query ", data("dummy_reads.fastq"),
                                                        "--ref-meta", ref_meta_path);

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
    auto const [number_of_bins, window_size, number_of_errors] = GetParam();

    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);
    size_t pattern_size = 50;
    float error_rate = (float) number_of_errors / (float) pattern_size;

    size_t query_seg_count = 60;
    std::filesystem::path ref_meta_path = "ref_meta.bin";
    std::filesystem::path index_path = ibf_path(number_of_bins, window_size);
    std::filesystem::path query_meta_path = "query_meta.bin";

    cli_test_result const split_ref = execute_app("valik", "split",
                                                        data("ref.fasta"),
                                                        "--out", ref_meta_path,
                                                        "--seg-count ", std::to_string(number_of_bins),
                                                        "--pattern 150", 
                                                        "--without-parameter-tuning");
    EXPECT_EQ(split_ref.exit_code, 0);

    cli_test_result const split_query = execute_app("valik", "split",
                                                        data("query.fasta"),
                                                        "--out", query_meta_path,
                                                        "--seg-count ", std::to_string(query_seg_count),
                                                        "--ref-meta", ref_meta_path,
                                                        "--pattern 0");
    EXPECT_EQ(split_query.exit_code, 0);
    valik::metadata reference(ref_meta_path);

    cli_test_result const build = execute_app("valik", "build",
                                                        "--window", std::to_string(window_size),
                                                        "--kmer 13",
                                                        "--size 32k",
                                                        "--ref-meta", ref_meta_path,
                                                        "--output ", index_path);
    EXPECT_EQ(build.exit_code, 0);

    cli_test_result const search = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--pattern ", std::to_string(pattern_size),
                                                        "--overlap ", std::to_string(pattern_size - 1),
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", index_path,
                                                        "--query ", data("query.fasta"),
                                                        "--ref-meta", ref_meta_path,
                                                        "--query-meta", query_meta_path, 
                                                        "--repeatPeriod 1",
                                                        "--repeatLength 10", 
                                                        "--numMatches 2");

    EXPECT_EQ(search.exit_code, 0);
    EXPECT_EQ(search.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 4 database sequences.\n"});
    EXPECT_EQ(search.err, std::string{});

    auto distributed = valik::read_stellar_output(search_result_path(number_of_bins, window_size, number_of_errors), reference, std::ios::binary);
    auto local = valik::read_stellar_output("search.gff", reference);
    
    compare_gff_out(distributed, local);
}

INSTANTIATE_TEST_SUITE_P(split_shared_memory_dream_suite,
                         dream_split_search,
                         testing::Combine(testing::Values(4, 16), testing::Values(13, 15), testing::Values(1)),
                         [] (testing::TestParamInfo<dream_split_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });
