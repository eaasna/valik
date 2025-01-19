#include <gtest/gtest.h>

#include <utilities/threshold/basics.hpp>

TEST(kmer_occurrences, small)
{
    uint8_t kmer_size = 4;
    uint64_t genome_len = 256 + kmer_size - 1;
    double inf_cont = 1.0;
    EXPECT_EQ(valik::expected_kmer_occurrences(genome_len, kmer_size, inf_cont), 1.0);
}

TEST(kmer_occurrences, human_genome)
{
    uint64_t human_genome_len = 3e9;
    uint8_t kmer_size = 16;
    double inf_cont = 1.0;
    EXPECT_LE(valik::expected_kmer_occurrences(human_genome_len, kmer_size, inf_cont), 1.0);
}

TEST(kmer_lemma, manual_test)
{
    size_t l = 100;
    uint8_t k = 10;

    uint8_t shared_kmer_count = l - k + 1; // 91
    for (uint8_t e{0}; e < 9; e++)
        EXPECT_EQ(valik::kmer_lemma_threshold(l, k, e), shared_kmer_count - k * e);

    uint8_t e = 10;
    EXPECT_EQ(valik::kmer_lemma_threshold(l, k, e), 0);

    l = 15;
    EXPECT_EQ(valik::kmer_lemma_threshold(l, k, e), 0);
}

TEST(combinations, edge_cases)
{
    uint8_t k = 3;
    size_t n = 3;
    EXPECT_EQ((valik::combinations(k, n)), 1);

    k = 4;
    EXPECT_EQ((valik::combinations(k, n)), 0);

    k = 0;
    n = 0;
    EXPECT_EQ((valik::combinations(k, n)), 1);
}

TEST(combinations, manual_test)
{
    uint8_t k = 3;
    size_t n = 4;
    EXPECT_EQ((valik::combinations(k, n)), 4);
    EXPECT_EQ((valik::combinations(k, n)), (valik::combinations(n - k, n)));

    n = 5;
    EXPECT_EQ((valik::combinations(k, n)), 10);
    EXPECT_EQ((valik::combinations(k, n)), (valik::combinations(n - k, n)));

    k = 13;
    n = 37;
    EXPECT_EQ((valik::combinations(k, n)), 3562467300); // check for numerical overflow
}
