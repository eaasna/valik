#include <gtest/gtest.h>

#include <filesystem>

#include <utilities/threshold/kmer_loss.hpp>
#include <utilities/threshold/io.hpp>
#include <utilities/threshold/search_pattern.hpp>
#include <utilities/threshold/param_set.hpp>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/core/debug_stream.hpp>

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data_path(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(kmer_loss, equal_after_serialization)
{
    std::filesystem::remove(valik::fn_filename());
    std::vector<valik::kmer_loss> precalc_attr_vec;
    valik::precalc_fn_confs(precalc_attr_vec);

    std::vector<valik::kmer_loss> deserialized_attr_vec;
    valik::read_fn_confs(deserialized_attr_vec);
    EXPECT_EQ(precalc_attr_vec.size(), deserialized_attr_vec.size());
    for (size_t i{0}; i < precalc_attr_vec.size(); i++)
    {
        EXPECT_EQ(precalc_attr_vec[i].k, deserialized_attr_vec[i].k);
        EXPECT_EQ(precalc_attr_vec[i].fn_conf_counts.size(), deserialized_attr_vec[i].fn_conf_counts.size());
        for (size_t j{0}; j < precalc_attr_vec[i].fn_conf_counts.size(); j++)
        {
            EXPECT_EQ(precalc_attr_vec[i].fn_conf_counts[j].size(), deserialized_attr_vec[i].fn_conf_counts[j].size());
            for (size_t k{0}; k < precalc_attr_vec[i].fn_conf_counts[j].size(); k++)
            {
                EXPECT_RANGE_EQ(precalc_attr_vec[i].fn_conf_counts[j][k], deserialized_attr_vec[i].fn_conf_counts[j][k]);
            }
        }
    }
}

TEST(kmer_loss, basic_checks)
{
    valik::param_space space;

    std::vector<valik::kmer_loss> attr_vec;
    if (!valik::read_fn_confs(attr_vec))
        valik::precalc_fn_confs(attr_vec);

    EXPECT_EQ(attr_vec.size(), std::get<1>(space.kmer_range) - std::get<0>(space.kmer_range) + 1);

    for (auto & attr : attr_vec)
    {
        EXPECT_EQ(attr.fn_conf_counts.size(), space.max_thresh);
        for (auto & thresh_table : attr.fn_conf_counts)
        {
            EXPECT_EQ(thresh_table.size(), space.max_errors + 1);
            for (auto & error_row : thresh_table)
            {
                EXPECT_EQ(error_row.size(), space.max_len + 1);
            }
        }
    }
}

// Test: case where there are no shared k-mers because the sequence length is less than the k-mer length
//       the conf count is equal to (l take e) because all possible error configuration destroy all shared k-mers
static void check_len_less_than_kmer_size(valik::kmer_loss const & attr)
{
    for (uint8_t t{1}; t <= attr.fn_conf_counts.size(); t++)
    {
        auto thresh_table = attr.fn_conf_counts[t - 1];
        for (size_t e{0}; e < thresh_table.size(); e++)
        {
            auto error_row = thresh_table[e];
            for (size_t l{0}; l < attr.k; l++)
            {
                if (error_row[l] != valik::combinations(e, l))
                {
                    std::cout << "l\t" << std::to_string(l) << '\n';
                    std::cout << "k\t" << std::to_string(attr.k) << '\n';
                    std::cout << "e\t" << std::to_string(e) << '\n';
                    std::cout << "t\t" << std::to_string(t) << '\n';
                }
                EXPECT_EQ(error_row[l], valik::combinations(e, l));
            }
        }
    }
}

// Test: case where the thresh parameter is less or equal to the k-mer lemma threshold
//       the conf count is 0 because no error configuration destroys more than the allowed number of k-mers
static void check_thresh_from_kmer_lemma(valik::kmer_loss const & attr)
{
    for (uint8_t t{1}; t <= attr.fn_conf_counts.size(); t++)
    {
        auto thresh_table = attr.fn_conf_counts[t - 1];
        for (uint8_t e{0}; e < thresh_table.size(); e++)
        {
            auto error_row = thresh_table[e];
            for (size_t l{attr.k}; l < error_row.size(); l++)
            {
                if (valik::kmer_lemma_threshold(l, attr.k, e) >= t)
                {
                    if (error_row[l] != 0)
                    {
                        std::cout << "l\t" << std::to_string(l) << '\n';
                        std::cout << "k\t" << std::to_string(attr.k) << '\n';
                        std::cout << "e\t" << std::to_string(e) << '\n';
                        std::cout << "t\t" << std::to_string(t) << '\n';
                    }
                    EXPECT_EQ(error_row[l], 0);
                }
            }
        }
    }
}

TEST(kmer_loss, edge_cases)
{
    std::vector<valik::kmer_loss> attr_vec;
    if (!valik::read_fn_confs(attr_vec))
        valik::precalc_fn_confs(attr_vec);

    for (auto & attr : attr_vec)
    {
        check_len_less_than_kmer_size(attr);
        check_thresh_from_kmer_lemma(attr);
    }
}

/*
    Let f(k, t, e, l) be the number of error configurations that destroy more than t k-mers. 
    The number of error configurations can be calculated in two ways.
*/
TEST(kmer_loss, exhaustive_comparison)
{
    std::vector<valik::kmer_loss> attr_vec;
    if (!valik::read_fn_confs(attr_vec))
        valik::precalc_fn_confs(attr_vec);

    valik::param_space space;
    for (uint8_t k{std::get<0>(space.kmer_range)}; k <= std::get<1>(space.kmer_range); k++)
    {   
        auto kmer_attr = attr_vec[k - 9];
        for (size_t t{1}; t <= space.max_thresh; t++)
        {
            for (size_t e{1}; e <= space.max_errors; e++)
            {
                for (size_t l = k + t - 1 + e; l <= space.max_len; l++)
                {
                    // f(k, t, e, l) is calculated based on f(k, t, e, l - 1)
                    uint64_t calc_conf_count = attr_vec[k - 9].fn_conf_counts[t - 1][e][l];

                    // f(k, t, e, l) is calculated based on a subvector sum in row (e - 1) 
                    // i.e f(k, t, e, l) = f(k, t, e - 1, l - 1) + f(k, t, e - 1, l - 2) ...     
                    uint64_t row_sum_conf_count{0};
                    for (size_t i{1}; i <= k; i++)
                        row_sum_conf_count += attr_vec[k - 9].fn_conf_counts[t - 1][e - 1][l - i];
        
                    for (size_t j{1}; j < t; j++)
                        row_sum_conf_count += attr_vec[k - 9].fn_conf_counts[t - 1 - j][e - 1][l - k - j];
    
                    EXPECT_EQ(row_sum_conf_count, calc_conf_count);
                }
            }
        }
    }
}

void try_fnr(uint8_t e, size_t l, uint8_t k, uint8_t t, double expected_fnr)
{
    valik::param_space space;
    valik::search_pattern pattern(e, l);
    valik::param_set param(k, t, space);

    std::vector<valik::kmer_loss> attr_vec;
    if (!valik::read_fn_confs(attr_vec))
        valik::precalc_fn_confs(attr_vec);

    EXPECT_EQ(attr_vec[k - std::get<0>(space.kmer_range)].fnr_for_param_set(pattern, param), expected_fnr); 
}

TEST(false_negative, try_kmer_lemma_thershold)
{
    try_fnr(0u, (size_t) 50, 16u, 35u, 0.0);
    try_fnr(1u, (size_t) 50, 16u, 19u, 0.0);
    try_fnr(2u, (size_t) 50, 16u, 3u, 0.0);
}

TEST(false_negative, try_threshold_above_kmer_lemma)
{
    uint8_t t{35};
    EXPECT_THROW({
        try
        {
            try_fnr(0u, (size_t) 15, 16u, t, 1.0);
        }
        catch( const std::runtime_error& e )
        {
            EXPECT_STREQ( ("Calculated configuration count table for t=[1, 5]. " 
                           "Can't find FNR for t=" + std::to_string(t)).c_str(), e.what() );
            throw;
        }
    }, std::runtime_error );
}

TEST(false_negative, try_all_destroyed)
{
    try_fnr(5u, (size_t) 10, 9u, 3u, 1.0);
    try_fnr(4u, (size_t) 10, 9u, 3u, 1.0);
    try_fnr(3u, (size_t) 10, 9u, 3u, 1.0);
}
