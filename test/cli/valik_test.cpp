#include <filesystem>
#include <fstream>
#include <ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include "cli_test.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// valik split preprocess clusters ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct valik_split_clusters : public valik_base {};

TEST_F(valik_split_clusters, split_metagenome_clusters)
{
    size_t genome_count{8};
    {
        std::ofstream one_per_bin{"single_seq_bin_paths.txt"};
        for (size_t i{0}; i < genome_count; i++)
        {
            std::string file_path = cli_test::data("bin_" + std::to_string(i) + ".fasta");
            one_per_bin << file_path << '\n';
        }
    
        std::ofstream two_per_bin{"multi_seq_bin_paths.txt"};
        for (size_t i{0}; i < genome_count; i = i + 2)
        {
            std::string file_path = cli_test::data("bin_" + std::to_string(i) + ".fasta");
            two_per_bin << file_path << '\t';
            file_path = cli_test::data("bin_" + std::to_string(i+1) + ".fasta");
            two_per_bin << file_path << '\n';
        }
    }

    cli_test_result const result_one_per_bin = execute_app("valik", "split",
                                                                    "--metagenome",
                                                                    "--split-index",
                                                                    "--out single_seq_meta.bin",
                                                                    "single_seq_bin_paths.txt");
    EXPECT_EQ(result_one_per_bin.exit_code, 0);
    EXPECT_EQ(result_one_per_bin.out, std::string{});
    EXPECT_EQ(result_one_per_bin.err, std::string{});

    auto one_per_bin_meta = valik::metadata("single_seq_meta.bin");
    EXPECT_EQ(one_per_bin_meta.seq_count, genome_count * 2);
    EXPECT_EQ(one_per_bin_meta.seg_count, genome_count);

    cli_test_result const result_two_per_bin = execute_app("valik", "split",
                                                                    "--metagenome",
                                                                    "--split-index",
                                                                    "--out multi_seq_meta.bin",
                                                                    "multi_seq_bin_paths.txt");
    EXPECT_EQ(result_two_per_bin.exit_code, 0);
    EXPECT_EQ(result_two_per_bin.out, std::string{});
    EXPECT_EQ(result_two_per_bin.err, std::string{});

    auto two_per_bin_meta = valik::metadata("multi_seq_meta.bin");
    EXPECT_EQ(two_per_bin_meta.seq_count, genome_count * 2);
    EXPECT_EQ(two_per_bin_meta.seg_count, genome_count / 2u);

    EXPECT_EQ(one_per_bin_meta.total_len, two_per_bin_meta.total_len);
    EXPECT_EQ(one_per_bin_meta.total_len, 8192*2);    // hard coded value from test/data/build/cli_test_input.sh
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// valik split index bins /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(valik_split_ref, split_ref)
{
    auto const [seg_count, overlap] = GetParam();

    cli_test_result const result = execute_app("valik", "split",
                                                        data("ref.fasta"),
                                                        "--split-index",
                                                        "--out reference_metadata.txt",
                                                        "--seg-count ", std::to_string(seg_count),
                                                        "--pattern ", std::to_string(overlap));

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});


    valik::metadata meta("reference_metadata.txt");
    if (meta.seg_count < 97)
        EXPECT_EQ(meta.seg_count, 64);
    else
        EXPECT_EQ(meta.seg_count, 128);
    EXPECT_GE(0.1f, meta.segment_length_cv()); // create segments of roughly equal length
}

INSTANTIATE_TEST_SUITE_P(split_ref_suite,
                         valik_split_ref,
                         testing::Combine(testing::Values(8, 63, 64, 65, 96, 97, 159), testing::Values(20, 21)),
                         [] (testing::TestParamInfo<valik_split_long::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_seg_count_" +
                                                std::to_string(std::get<1>(info.param)) + "_overlap";
                             return name;
                         });


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// valik split equal length ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(valik_split_various, split_various_lengths)
{
    auto const [seg_count, overlap] = GetParam();

    cli_test_result const result = execute_app("valik", "split",
                                                        data("various_chromosome_lengths.fasta"),
                                                        "--out query_metadata.txt",
                                                        "--seg-count ", std::to_string(seg_count),
                                                        "--pattern ", std::to_string(overlap), 
                                                        "--ref-meta ", segment_metadata_path(150, 4));

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"Sequence: chr5 is too short and will be skipped.\n"}); 
    auto expected_segments = valik::metadata(segment_metadata_path(overlap, seg_count));
    std::string expected_segment_str = expected_segments.to_string();
    auto actual_segments = valik::metadata("query_metadata.txt");
    std::string actual_segment_str = actual_segments.to_string();
    EXPECT_TRUE(expected_segment_str == actual_segment_str);
}


INSTANTIATE_TEST_SUITE_P(split_suite,
                         valik_split_various,
                         testing::Combine(testing::Values(4, 16), testing::Values(20)),
                         [] (testing::TestParamInfo<valik_split_various::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_seg_count_" +
                                                std::to_string(std::get<1>(info.param)) + "_overlap";
                             return name;
                         });

TEST_P(valik_split_short, split_many_short)
{
    auto const [seg_count, overlap] = GetParam();

    cli_test_result const result = execute_app("valik", "split",
                                                        data("query.fasta"),
                                                        "--out query_metadata.txt",
                                                        "--ref-meta ", segment_metadata_path(150, 4), 
                                                        "--seg-count ", std::to_string(seg_count),
                                                        "--pattern ", std::to_string(overlap));

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});

    size_t read_count = 31;
    if (seg_count % read_count == 0)
    {
        EXPECT_EQ(result.err, std::string{});
    }
    // it is not advantageous to split 30 sequences of equal length into e.g 39 segments
    else
    {
        const std::string expected = "WARNING: Database was split into " +
                                     std::to_string((size_t) std::round((float) seg_count / read_count) * read_count) +
                                     " instead of " + std::to_string(seg_count) + " segments.\n";
        EXPECT_EQ(result.err, expected);
    }

    valik::metadata meta("query_metadata.txt");
    EXPECT_GE(0.1f, meta.segment_length_cv());  // create segments of roughly equal length
}

INSTANTIATE_TEST_SUITE_P(split_many_short_suite,
                         valik_split_short,
                         testing::Combine(testing::Values(31, 39, 41, 49, 51, 55, 60, 61, 71), testing::Values(20)),
                         [] (testing::TestParamInfo<valik_split_short::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_seg_count_" +
                                                std::to_string(std::get<1>(info.param)) + "_overlap";
                             return name;
                         });

TEST_P(valik_split_long, split_few_long)
{
    auto const [seg_count, overlap] = GetParam();

    cli_test_result const result = execute_app("valik", "split",
                                                        data("ref.fasta"),
                                                        "--out reference_metadata.txt",
                                                        "--ref-meta", segment_metadata_path(150, 4),
                                                        "--seg-count ", std::to_string(seg_count),
                                                        "--pattern ", std::to_string(overlap));

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});

    valik::metadata meta("reference_metadata.txt");
     
    if (seg_count == meta.seq_count) // one-to-one pairing of sequences and segments
    {
        EXPECT_EQ(result.err, std::string{"WARNING: Database was split into 5 instead of 4 segments.\n"});
    }
    else
    {
        EXPECT_GE(0.2f, meta.segment_length_cv());  // create segments of roughly equal length
        EXPECT_EQ(result.err, std::string{});
    }
}

INSTANTIATE_TEST_SUITE_P(split_few_long_suite,
                         valik_split_long,
                         testing::Combine(testing::Values(4, 12, 19), testing::Values(20)),
                         [] (testing::TestParamInfo<valik_split_long::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_seg_count_" +
                                                std::to_string(std::get<1>(info.param)) + "_overlap";
                             return name;
                         });

struct split_options : public valik_base {};

TEST_F(split_options, too_few_segments)
{
    size_t n = 30;
    size_t o = 0;
    cli_test_result const result = execute_app("valik", "split", data("query.fasta"), 
                                                        "--ref-meta", segment_metadata_path(150, 4),
                                                        "--seg-count", std::to_string(n), "--pattern", std::to_string(o),
                                                        "--out", "meta.txt");
    std::string const expected
    {
        "[Error] Can not split 31 sequences into " + std::to_string(n) + " segments.\n"
    };

    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(split_options, overlap_too_large)
{
    size_t n = 30;
    size_t o = 2000;
    cli_test_result const result = execute_app("valik", "split", data("query.fasta"), 
                                                        "--ref-meta", segment_metadata_path(150, 4),
                                                        "--seg-count", std::to_string(n), "--pattern", std::to_string(o),
                                                        "--out", "meta.txt");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    std::string err_first_half = result.err.substr(0, result.err.find("length") + 6);
    std::string err_second_half = result.err.substr(result.err.find("bp") + 2, result.err.size() - 1);
    EXPECT_EQ(err_first_half, std::string{"[Error] Segments of length"});
    std::string expected = " can not overlap by " + std::to_string(o) + "bp.\nDecrease the overlap or the number of segments.\n";
    EXPECT_EQ(err_second_half, expected);
}

TEST_F(split_options, too_many_segments)
{
    size_t n = 300;
    size_t o = 20;
    cli_test_result const result = execute_app("valik", "split", data("query.fasta"), 
                                                        "--ref-meta", segment_metadata_path(150, 4),
                                                        "--seg-count", std::to_string(n), "--pattern", std::to_string(o),
                                                        "--out", "meta.txt");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    std::string err_first_half = result.err.substr(0, result.err.find("length") + 6);
    std::string err_second_half = result.err.substr(result.err.find("bp") + 2, result.err.size() - 1);
    EXPECT_EQ(err_first_half, std::string{"[Error] Segments of length"});
    std::string expected = " can not overlap by " + std::to_string(o) + "bp.\nDecrease the overlap or the number of segments.\n";
    EXPECT_EQ(err_second_half, expected);
}

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
    
    std::filesystem::path ref_meta_path = "bin_meta.bin";
    cli_test_result const split_ref = execute_app("valik", "split",
                                                           "bin_paths.txt",
                                                           "--out", ref_meta_path,
                                                           "--metagenome");
    EXPECT_EQ(split_ref.exit_code, 0);

    cli_test_result const result = execute_app("valik", "build",
                                                         "--kmer 19",
                                                         "--window ", std::to_string(window_size),
                                                         "--size 32k",
                                                         "--threads ", run_parallel ? "2" : "1",
                                                         "--output index.ibf", 
                                                         "--ref-meta ", ref_meta_path);
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

    std::filesystem::path ref_meta_path = "ref_meta.bin";
    cli_test_result const split_ref = execute_app("valik", "split",
                                                           data("single_reference.fasta"),
                                                           "--out", ref_meta_path,
                                                           "--pattern 150", 
                                                           "--seg-count ", std::to_string(number_of_bins), 
                                                           "--without-parameter-tuning");
    EXPECT_EQ(split_ref.exit_code, 0);
    EXPECT_EQ(split_ref.out, std::string{});

    cli_test_result const result = execute_app("valik", "build",
                                                        "--kmer 13",
                                                        "--window ", std::to_string(window_size),
                                                        "--size 32k",
                                                        "--output index.ibf",
                                                        "--ref-meta", ref_meta_path);
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
    auto const [number_of_bins, window_size, number_of_errors, pattern_size, overlap] = GetParam();
    float error_rate = (float) number_of_errors / (float) pattern_size;

    setup_tmp_dir();
    setenv("VALIK_STELLAR", "echo", true);
    setenv("VALIK_MERGE", "echo", true);

    if (window_size == 23 && number_of_errors == 0)
        GTEST_SKIP() << "Needs dynamic threshold correction";

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--distribute",
                                                        "--pattern", std::to_string(pattern_size),
                                                        "--overlap", std::to_string(overlap),
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", ibf_path(number_of_bins, window_size),
                                                        "--query ", data("query.fq"),
                                                        "--threads 1",
                                                        "--tau 0.75",
                                                        "--p_max 0.75");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    auto expected = string_list_from_file(search_result_path(number_of_bins, window_size, number_of_errors,
                                                             pattern_size, overlap), std::ios::binary);
    auto actual = string_list_from_file("search.gff.out");

    EXPECT_EQ(expected, actual);
}

INSTANTIATE_TEST_SUITE_P(cluster_search_suite,
                         valik_search_clusters,
                         testing::Combine(testing::Values(8), testing::Values(19, 23), testing::Values(0, 1),
                 testing::Values(50, 100), testing::Values(1, 40)),
                         [] (testing::TestParamInfo<valik_search_clusters::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error_" +
                                                std::to_string(std::get<3>(info.param)) + "_pattern_" +
                                                std::to_string(std::get<4>(info.param)) + "_overlap";
                             return name;
                         });

TEST_P(valik_search_segments, search)
{
    auto const [segment_overlap, number_of_bins, window_size, number_of_errors, pattern_size, overlap] = GetParam();
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
                                                        "--overlap", std::to_string(overlap),
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", ibf_path(segment_overlap, number_of_bins, window_size),
                                                        "--query ", data("single_query.fasta"),
                                                        "--tau 0.75",
                                                        "--threads 1",
                                                        "--ref-meta", segment_metadata_path(segment_overlap, number_of_bins),
                                                        "--p_max 0.25");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    auto expected = string_list_from_file(search_result_path(segment_overlap, number_of_bins, window_size, 
                                                             number_of_errors), std::ios::binary);
    auto actual = string_list_from_file("search.gff.out");

    EXPECT_EQ(expected, actual);
}

INSTANTIATE_TEST_SUITE_P(segment_search_suite,
                         valik_search_segments,
                         testing::Combine(testing::Values(150), testing::Values(4, 16), testing::Values(15, 13), testing::Values(1),
                         testing::Values(50), testing::Values(49)),
                         [] (testing::TestParamInfo<valik_search_segments::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<1>(info.param)) + "_bins_" +
                                                std::to_string(std::get<2>(info.param)) + "_window_" +
                                                std::to_string(std::get<3>(info.param)) + "_error_" +
                                                std::to_string(std::get<4>(info.param)) + "_pattern_" +
                                                std::to_string(std::get<5>(info.param)) + "_overlap";
                             return name;
                         });
