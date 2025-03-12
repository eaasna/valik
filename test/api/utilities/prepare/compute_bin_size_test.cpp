#include <gtest/gtest.h>

#include <utilities/prepare/compute_bin_size.hpp>
#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>
#include <valik/split/write_seg_sequences.hpp>

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(max_bin_count, split_db_no_count_cutoff)
{
    valik::build_arguments arguments{};
    arguments.bin_path = std::vector<std::string>{};
    size_t bin_count = 8;
    for (size_t b{0}; b < bin_count; b++)
        arguments.bin_path.emplace_back(data("ref_bin_" + std::to_string(b) + ".fasta"));

    std::vector<seqan3::shape> shapes{seqan3::shape{seqan3::ungapped{8}}, seqan3::shape{seqan3::bin_literal{0b1001}}};
    for (auto shape : shapes)
    {
        arguments.shape = shape;
        arguments.window_size = shape.size() + 2;
        auto sequence_files_max_count = raptor::detail::kmer_count_from_sequence_files(arguments);
        std::vector<std::string> header_files{};
        for (size_t b{0}; b < bin_count; b++)
            header_files.emplace_back(data("s" + shape.to_string() + "_ref." + std::to_string(b) + ".header"));

        auto minimiser_files_max_count = raptor::detail::kmer_count_from_minimiser_files(header_files, arguments.threads);
    
        EXPECT_EQ(sequence_files_max_count, minimiser_files_max_count);    
    }
}

TEST(max_bin_count, segment_db_no_count_cutoff)
{
    size_t bin_count = 8;
    valik::build_arguments arguments{};
    arguments.bin_path = std::vector<std::string>{data("ref.fasta")};
    
    std::vector<seqan3::shape> shapes{seqan3::shape{seqan3::ungapped{8}}, seqan3::shape{seqan3::bin_literal{0b1001}}};
    for (auto shape : shapes)
    {
        arguments.ref_meta_path = data("s" + shape.to_string() + ".bin");
        arguments.shape = shape;
        arguments.window_size = shape.size() + 2;
        auto sequence_files_max_count = raptor::detail::kmer_count_from_sequence_files(arguments);
        std::vector<std::string> header_files{};
        for (size_t b{0}; b < bin_count; b++)
            header_files.emplace_back(data("s" + shape.to_string() + "_ref." + std::to_string(b) + ".header"));

        auto minimiser_files_max_count = raptor::detail::kmer_count_from_minimiser_files(header_files, arguments.threads);
    
        EXPECT_EQ(sequence_files_max_count, minimiser_files_max_count);    
    }
}

/*
TEST(compute_bin_minimisers, segment_db_no_count_cutoff)
{
    valik::split_arguments split_arguments{};
    split_arguments.pattern_size = 50;
    split_arguments.bin_path = std::vector<std::string>{data("ref.fasta")};
    split_arguments.seg_count_in = 8;    
    split_arguments.write_out = true;


    valik::build_arguments build_arguments{};
    std::vector<size_t> bin_ids{0, 20, 40, 60};
    build_arguments.bin_path = std::vector<std::string>{data("ref.fasta")};
    std::vector<seqan3::shape> shapes{seqan3::shape{seqan3::ungapped{8}}, seqan3::shape{seqan3::bin_literal{0b1001}}};
    for (auto shape : shapes)
    {
        split_arguments.shape = shape;
        split_arguments.kmer_size = shape.size();
        split_arguments.shape_str = shape.to_string();
        split_arguments.meta_out = "s" + shape.to_string() + ".bin";
        split_arguments.shape_weight = shape.count();
        auto meta = valik::metadata{split_arguments};
        meta.save(split_arguments.meta_out);
        valik::write_reference_segments(meta, data("ref.fasta"));
        
        build_arguments.ref_meta_path = data(split_arguments.meta_out);
        build_arguments.shape = shape;
        build_arguments.window_size = shape.size() + 2;
        raptor::compute_minimiser(build_arguments);
        auto sequence_files_max_count = raptor::detail::kmer_count_from_sequence_files(build_arguments);        
    }
}
*/