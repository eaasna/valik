#include <gtest/gtest.h>

#include <valik/split/metadata.hpp>

// Generate the full path of a test input file that is provided in the data directory.
static std::filesystem::path data_path(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

static std::filesystem::path segment_metadata_path(size_t const overlap, size_t const bins) noexcept
{
    std::string name{};
    name += std::to_string(overlap);
    name += "overlap";
    name += std::to_string(bins);
    name += "bins.bin";
    return data_path(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// valik split preprocess clusters ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(split_options, split_metagenome_clusters)
{
    valik::build_arguments arguments{};
    arguments.metagenome = true;

    size_t genome_count{8};
    {
        for (size_t i{0}; i < genome_count; i++)
        {
            std::string file_path = data_path("bin_" + std::to_string(i) + ".fasta");
            arguments.bin_path.emplace_back(file_path);
        }
    }

    try
    {
        valik::metadata meta(arguments);
        EXPECT_EQ(meta.seq_count, genome_count * 2);
        EXPECT_EQ(meta.seg_count, genome_count);
        EXPECT_EQ(meta.total_len, 8192*2);    // hard coded value from test/data/build/cli_test_input.sh
    }
    catch( const std::runtime_error& e )
    {
        EXPECT_TRUE( false );
        throw;
    }    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// valik split index bins /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(split_options, split_ref)
{
    valik::build_arguments arguments{};
    arguments.error_rate = 0.01;
    arguments.bin_path.emplace_back(data_path("ref.fasta"));

    try
    {
        for (auto b : std::vector<uint32_t>{8, 63, 64, 65, 96, 97, 159})
        {
            arguments.bins = b;
            for (auto p : std::vector<uint8_t>{20, 21})
            {
                arguments.pattern_size = p;
                valik::metadata meta(arguments);
                if (b < 97UL)
                    EXPECT_EQ(meta.seg_count, 64UL);
                else
                    EXPECT_EQ(meta.seg_count, 128UL);
                EXPECT_GE(0.1f, meta.segment_length_cv()); // create segments of roughly equal length
            }
        }
    }
    catch( const std::runtime_error& e )
    {
        EXPECT_TRUE( false );
        throw;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// valik split equal length ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(split_options, split_various_length)
{
    valik::build_arguments arguments{};
    arguments.error_rate = 0.01;
    arguments.only_split = true;
    arguments.pattern_size = 20;
    arguments.bin_path.emplace_back(data_path("various_chromosome_lengths.fasta"));

    try
    {
        for (arguments.bins : std::vector<uint32_t>{4, 16})
        {
            valik::metadata meta(arguments);
            auto expected_segments = valik::metadata(segment_metadata_path(arguments.pattern_size, arguments.bins));
            std::string expected_segment_str = expected_segments.to_string();
            auto actual_segments = valik::metadata("various_chromosome_lengths.bin");
            std::string actual_segment_str = actual_segments.to_string();
            EXPECT_TRUE(expected_segment_str == actual_segment_str);            
        }
    }
    catch( const std::runtime_error& e )
    {
        EXPECT_TRUE( false );
        throw;
    }
}

TEST(split_options, split_few_long)
{
    valik::build_arguments arguments{};
    arguments.error_rate = 0.01;
    arguments.only_split = true;
    arguments.pattern_size = 20;
    arguments.bin_path.emplace_back(data_path("various_chromosome_lengths.fasta"));

    try
    {
        for (arguments.bins : std::vector<uint32_t>{4, 12, 19})
        {
            valik::metadata meta(arguments);

            if (seg_count > meta.seq_count) // one-to-many pairing of sequences and segments
            {
                EXPECT_GE(0.2f, meta.segment_length_cv());  // create segments of roughly equal length
            }
        }
    }
    catch( const std::runtime_error& e )
    {
        EXPECT_TRUE( false );
        throw;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// valik split error handling ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(split_options, too_few_segments)
{
    valik::build_arguments arguments{};
    arguments.bin_path.emplace_back(data_path("query.fasta"));
    arguments.seg_count = 30;
    arguments.pattern_size = 10;

    EXPECT_THROW({
        try
        {
            valik::metadata meta(arguments);
        }
        catch( const std::runtime_error& e )
        {
            EXPECT_STREQ( ("Can not split 31 sequences into " + std::to_string(arguments.seg_count) + " segments.").c_str(), e.what() );
            throw;
        }
    }, std::runtime_error );
}

TEST(split_options, overlap_too_large)
{
    valik::build_arguments arguments{};
    arguments.bin_path.emplace_back(data_path("query.fasta"));
    arguments.seg_count = 30;
    arguments.pattern_size = 1000;

    EXPECT_THROW({
        try
        {
            valik::metadata meta(arguments);
        }
        catch( const std::runtime_error& e )
        {
            std::string err = e.what();
            std::string err_first_half = err.substr(0, err.find("length") + 6);
            std::string err_second_half = err.substr(err.find("bp") + 2, err.size() - 1);
            EXPECT_EQ(err_first_half, "Segments of length");
            std::string expected = " can not overlap by " + std::to_string(arguments.pattern_size) + "bp.\nDecrease the overlap or the number of segments.";
            EXPECT_EQ(err_second_half, expected);

            throw;
        }
    }, std::runtime_error );
}

TEST(split_options, too_many_segment)
{
    valik::build_arguments arguments{};
    arguments.bin_path.emplace_back(data_path("query.fasta"));
    arguments.seg_count = 300;
    arguments.pattern_size = 20;

    EXPECT_THROW({
        try
        {
            valik::metadata meta(arguments);
        }
        catch( const std::runtime_error& e )
        {
            std::string err = e.what();
            std::string err_first_half = err.substr(0, err.find("length") + 6);
            std::string err_second_half = err.substr(err.find("bp") + 2, err.size() - 1);
            EXPECT_EQ(err_first_half, std::string{"Segments of length"});
            std::string expected = " can not overlap by " + std::to_string(arguments.pattern_size) + "bp.\nDecrease the overlap or the number of segments.";
            EXPECT_EQ(err_second_half, expected);
        
            throw;
        }
    }, std::runtime_error );
}
