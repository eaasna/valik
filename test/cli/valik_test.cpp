#include <filesystem>
#include <fstream>
#include <ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include "cli_test.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// valik build clusters /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(valik_build_clusters, build_from_clusters)
{
    auto const [number_of_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_bins >= 32;

    {
        std::ofstream file{"bin_paths.txt"};
        for (size_t i{0}; i < number_of_bins; i++)
        {
            std::string file_path = cli_test::data("bin_" + std::to_string(i) + ".fasta");
            file << file_path << '\n';
        }
        file << '\n';
    }

    std::cerr << "bin_paths\n";
    cli_test_result const result = execute_app("valik", "build",
                                                         "bin_paths.txt",
                                                         "--metagenome",
                                                         "--kmer 19",
                                                         "--window ", std::to_string(window_size),
                                                         "--size 32k",
                                                         "--threads ", run_parallel ? "2" : "1",
                                                         "--output index.ibf", 
                                                         "--without-parameter-tuning");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    compare_index(ibf_path(number_of_bins, window_size), "index.ibf");
}


INSTANTIATE_TEST_SUITE_P(cluster_build_suite,
                         valik_build_clusters,
                         testing::Combine(testing::Values(8), testing::Values(19, 23), testing::Values(true, false)),
                         [] (testing::TestParamInfo<valik_build_clusters::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                (std::get<2>(info.param) ? "parallel" : "serial");
                             return name;
                         });

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// valik build segments /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(valik_build_segments, build_from_segments)
{
    auto const [overlap, number_of_bins, window_size] = GetParam();

    //!TODO: faults because bin_paths is incorrect?

    cli_test_result const result = execute_app("valik", "build",
                                                        data("single_reference.fasta"),
                                                        "--kmer 13",
                                                        "--pattern 150", 
                                                        "--seg-count ", std::to_string(number_of_bins), 
                                                        "--window ", std::to_string(window_size),
                                                        "--size 32k",
                                                        "--output index.ibf", 
                                                        "--without-parameter-tuning");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    compare_index(ibf_path(overlap, number_of_bins, window_size), "index.ibf");
}


INSTANTIATE_TEST_SUITE_P(segment_build_suite,
                         valik_build_segments,
                         testing::Combine(testing::Values(150), testing::Values(4, 16), testing::Values(15, 13)),
                         [] (testing::TestParamInfo<valik_build_segments::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_overlap_" +
                                                std::to_string(std::get<1>(info.param)) + "_bins_" +
                                                std::to_string(std::get<2>(info.param)) + "_window";
                             return name;
                         });

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// valik distributed search ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(valik_search_clusters, search)
{
    auto const [number_of_bins, window_size, number_of_errors, pattern_size, query_every] = GetParam();
    float error_rate = (float) number_of_errors / (float) pattern_size;

    setup_tmp_dir();
    setenv("VALIK_STELLAR", "echo", true);
    setenv("VALIK_MERGE", "echo", true);

    if (window_size == 23 && number_of_errors == 0)
        GTEST_SKIP() << "Needs dynamic threshold correction";

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--pattern", std::to_string(pattern_size),
                                                        "--query-every", std::to_string(query_every),
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", ibf_path(number_of_bins, window_size),
                                                        "--query ", data("query.fq"),
                                                        "--threads 1", "--very-verbose",
                                                        "--cart-max-capacity 3",
                                                        "--max-queued-carts 10",
                                                        "--without-parameter-tuning");
    EXPECT_EQ(result.exit_code, 0);

    auto expected = string_list_from_file(search_result_path(number_of_bins, window_size, number_of_errors,
                                                             pattern_size, query_every), std::ios::binary);
    auto actual = string_list_from_file("search.gff.out");

    EXPECT_EQ(expected, actual);
}

INSTANTIATE_TEST_SUITE_P(cluster_search_suite,
                         valik_search_clusters,
                         testing::Combine(testing::Values(8), testing::Values(19, 23), testing::Values(0, 1),
                 testing::Values(50, 100), testing::Values(1, 50)),
                         [] (testing::TestParamInfo<valik_search_clusters::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error_" +
                                                std::to_string(std::get<3>(info.param)) + "_pattern_" +
                                                std::to_string(std::get<4>(info.param)) + "_query_every";
                             return name;
                         });

TEST_P(valik_search_segments, search)
{
    auto const [segment_overlap, number_of_bins, window_size, number_of_errors, pattern_size, query_every] = GetParam();
    float error_rate = (float) number_of_errors / (float) pattern_size;

    setup_tmp_dir();
    setenv("VALIK_STELLAR", "echo", true);
    setenv("VALIK_MERGE", "echo", true);

    // create dummy file for test case because VALIK_STELLAR == echo
    std::filesystem::path file_path{"search.gff.preliminary"};
    std::ofstream file_handle(file_path);
    file_handle << "\n";

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--distribute",
                                                        "--pattern", std::to_string(pattern_size),
                                                        "--query-every", std::to_string(query_every),
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", ibf_path(segment_overlap, number_of_bins, window_size),
                                                        "--query ", data("single_query.fasta"),
                                                        "--threads 1", "--very-verbose",
                                                        "--cart-max-capacity 3",
                                                        "--max-queued-carts 10",
                                                        "--without-parameter-tuning");
    EXPECT_EQ(result.exit_code, 0);

    auto expected = string_list_from_file(search_result_path(segment_overlap, number_of_bins, window_size, 
                                                             number_of_errors), std::ios::binary);
    auto actual = string_list_from_file("search.gff.out");

    EXPECT_EQ(expected, actual);
}

INSTANTIATE_TEST_SUITE_P(segment_search_suite,
                         valik_search_segments,
                         testing::Combine(testing::Values(150), testing::Values(4, 16), testing::Values(15, 13), testing::Values(1),
                         testing::Values(50), testing::Values(1)),
                         [] (testing::TestParamInfo<valik_search_segments::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<1>(info.param)) + "_bins_" +
                                                std::to_string(std::get<2>(info.param)) + "_window_" +
                                                std::to_string(std::get<3>(info.param)) + "_error_" +
                                                std::to_string(std::get<4>(info.param)) + "_pattern_" +
                                                std::to_string(std::get<5>(info.param)) + "_query_every";
                             return name;
                         });
