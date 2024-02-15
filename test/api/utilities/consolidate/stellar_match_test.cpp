#include <gtest/gtest.h>

#include <utilities/consolidate/stellar_match.hpp>

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data_path(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(print_match, no_evalue)
{
    std::vector<std::string> match_vec{"chr1", "Stellar", "eps-matches", "10", "108", "0.987", "+", ".", "1;seq2Range=1280,1378;cigar=97M1D2M;mutations=14A,45G,58T,92C"};
    std::string expected = "chr1\tStellar\teps-matches\t10\t108\t0.987\t+\t.\t1;seq2Range=1280,1378;cigar=97M1D2M;mutations=14A,45G,58T,92C\n";
    valik::metadata meta(data_path("8bins50overlap_reference_metadata.bin"));
    auto match = valik::stellar_match(match_vec, meta);
    std::string match_str = match.to_string();
    EXPECT_EQ(match_str, expected);
}

TEST(print_match, full)
{
    std::vector<std::string> match_vec{"chr1", "Stellar", "eps-matches", "10", "108", "0.987", "+", ".", "1;seq2Range=1280,1378;eValue=4.05784e-73;cigar=97M1D2M;mutations=14A,45G,58T,92C"};
    std::string expected = "chr1\tStellar\teps-matches\t10\t108\t0.987\t+\t.\t1;seq2Range=1280,1378;eValue=4.05784e-73;cigar=97M1D2M;mutations=14A,45G,58T,92C\n";
    valik::metadata meta(data_path("8bins50overlap_reference_metadata.bin"));
    auto match = valik::stellar_match(match_vec, meta);
    std::string match_str = match.to_string();
    EXPECT_EQ(match_str, expected);
}

TEST(create_match, malformed)
{
    std::vector<std::string> match_vec{"chr1", "Stellar", "eps-matches", "10", "108", "0.987", "+", ".", "1;seq2Range=1280,1378;eValue=4.05784e-73"};
    std::string attributes = "1;seq2Range=1280,1378;eValue=4.05784e-73";
    valik::metadata meta(data_path("8bins50overlap_reference_metadata.bin"));
    EXPECT_THROW({
        try
        {
            auto match = valik::stellar_match(match_vec, meta);
        }
        catch( const std::runtime_error& e )
        {
            // and this tests that it has the correct message
            EXPECT_EQ( "Malformed GFF record:\n" + attributes, e.what() );
            throw;
        }
    }, std::runtime_error );
}
