#include <gtest/gtest.h>

#include <utilities/prepare/compute_bin_size.hpp>
#include <valik/shared.hpp>

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(from_sequence_files, max_count)
{
    valik::build_arguments arguments{};
    arguments.bin_path = std::vector<std::string>{data("ref_bin_0.fasta"), data("ref_bin_20.fasta"), data("ref_bin_40.fasta"), data("ref_bin_60.fasta")};

    std::vector<seqan3::shape> shapes{seqan3::shape{seqan3::ungapped{8}}, seqan3::shape{seqan3::bin_literal{0b1001}}};
    for (auto shape : shapes)
    {
        arguments.shape = shape;
        arguments.window_size = shape.size() + 2;
        auto sequence_files_max_count = raptor::detail::kmer_count_from_sequence_files(arguments);
        auto minimiser_files_max_count = raptor::detail::kmer_count_from_minimiser_files(std::vector<std::string>{data("s" + shape.to_string() + "_ref.0.header"), 
                                                                                                                  data("s" + shape.to_string() + "_ref.20.header"), 
                                                                                                                  data("s" + shape.to_string() + "_ref.40.header"), 
                                                                                                                  data("s" + shape.to_string() + "_ref.60.header")},
                                                                                        arguments.threads);
    
        EXPECT_EQ(sequence_files_max_count, minimiser_files_max_count);    
    }
}
