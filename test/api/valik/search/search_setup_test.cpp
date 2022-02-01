#include <gtest/gtest.h>

#include <valik/search/search_setup.hpp>
#include <valik/shared.hpp>
#include <valik/search/compute_simple_model.hpp>

TEST(pattern_begin_positions, read_length_and_pattern_size_are_equal)
{
    // edge case where read_len = pattern_size
    size_t const read_len = 150u;
    uint64_t const pattern_size = read_len;
    uint64_t const overlap = 30u;

    std::vector<size_t> begin_positions{};
    std::vector<size_t> expected{0u}; // pattern interval [0, 150]

    valik::pattern_begin_positions(read_len, pattern_size, overlap, [&](size_t const begin)
    {
        begin_positions.push_back(begin);
    });

    EXPECT_EQ(begin_positions, expected); // seen all positions
}

TEST(pattern_begin_positions, overlaps_are_evenly)
{
    // special case where read_len - pattern_size is exactly the last possible begin position
    size_t const read_len = 150u;
    uint64_t const pattern_size = 50u;
    uint64_t const overlap = 30u;

    std::vector<size_t> begin_positions{};
    std::vector<size_t> expected = {
        0u, // pattern interval [0, 50]
        0u + pattern_size - overlap /* 20u */, // pattern interval [20, 70]
        20u + pattern_size - overlap /* 40u */, // pattern interval [40, 90]
        60u, // pattern interval [60, 110]
        80u, // pattern interval [80, 130]
        // pattern interval [100, 150]
        100u // ends here, because read_len - pattern_size = 100u
    };

    valik::pattern_begin_positions(read_len, pattern_size, overlap, [&](size_t const begin)
    {
        begin_positions.push_back(begin);
    });

    EXPECT_EQ(begin_positions, expected); // seen all positions
}

TEST(pattern_begin_positions, extra_overlap)
{
    // normal case where read_len - pattern_size has room for an additional element
    size_t const read_len = 150;
    uint64_t const pattern_size = 50;
    uint64_t const overlap = 20;

    std::vector<size_t> begin_positions{};
    std::vector<size_t> expected = {
        0u, // pattern interval [0, 50]
        0u + pattern_size - overlap /* 30u */, // pattern interval [30, 80]
        30u + pattern_size - overlap /* 60u */, // pattern interval [60, 110]
        90u, // pattern interval [90, 140]
        // pattern interval [100, 150]
        100u // ends here, because read_len - pattern_size = 110u
    };

    valik::pattern_begin_positions(read_len, pattern_size, overlap, [&](size_t const begin)
    {
        begin_positions.push_back(begin);
    });

    EXPECT_EQ(begin_positions, expected); // seen all positions
}

// ====================================================================================
// seq = CGCAAAACGCGGC      minimisers = [AAAA, AAAC, AACG]
// 	w = 8
// 	k = 4
//
// minimiser	             span_begin    minimiser belongs to windows starting at
// AAAA		CGCAAAACGCG          0                      0, 1, 2, 3
// AAAC		    AAACGCGG         4                          4
// AACG		     AACGCGGC        5                          5
//
// ====================================================================================

TEST(make_pattern_bounds, first_pattern_of_query)
{
    // normal case where pattern starts from beginning of query
    valik::search_arguments arguments{};
    arguments.pattern_size = 12;
    arguments.window_size = 8;
    arguments.kmer_size = 4;
    arguments.errors = 1;

    std::vector<size_t> const & window_span_begin{0, 4, 5};

    size_t const pattern_begin = 0;
    valik::pattern_bounds expected{};


    // |CGCAAAACGCGG|C
    //
    // 0: |[CGCAAAAC]GCGG|C     consider
    // 1: |C[GCAAAACG]CGG|C     consider
    // 2: |CG[CAAAACGC]GG|C     consider
    // 3: |CGC[AAAACGCG]G|C     consider
    // 4: |CGCA[AAACGCGG]|C     consider
    // 5: |CGCAA[AACGCGG |C]     don't consider
    //
    // minimisers[0, 2)
    expected.begin_position = 0;
    expected.end_position = 2;

    valik::threshold threshold_data = valik::make_threshold_data(arguments);
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, threshold_data);

    EXPECT_EQ(bounds.begin_position, expected.begin_position);
    EXPECT_EQ(bounds.end_position, expected.end_position);
}

TEST(make_pattern_bounds, same_minimiser_consecutive_windows_begin)
{
    // special case where pattern starts after the first of multiple consecutive minimisers
    // and pattern ends at the end of the query
    valik::search_arguments arguments{};
    arguments.pattern_size = 12;
    arguments.window_size = 8;
    arguments.kmer_size = 4;
    arguments.errors = 1;


    std::vector<size_t> const & window_span_begin{0, 4, 5};

    size_t const pattern_begin = 1;
    valik::pattern_bounds expected{};

    // C|GCAAAACGCGGC|
    //
    // 0: [C|GCAAAAC]GCGGC |     don't consider
    // 1:  C|[GCAAAACG]CGGC|     consider
    // 2:  C|G[CAAAACGC]GGC|     consider
    // 3:  C|GC[AAAACGCG]GC|     consider
    // 4:  C|GCA[AAACGCGG]C|     consider
    // 5:  C|GCAA[AACGCGGC]|     consider
    //
    // minimisers[0, 3)
    expected.begin_position = 0;
    expected.end_position = 3;

    valik::threshold threshold_data = valik::make_threshold_data(arguments);
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, threshold_data);

    EXPECT_EQ(bounds.begin_position, expected.begin_position);
    EXPECT_EQ(bounds.end_position, expected.end_position);
}

TEST(make_pattern_bounds, pattern_equals_window)
{
    // special case where pattern_size == window_size
    valik::search_arguments arguments{};
    arguments.pattern_size = 8;
    arguments.window_size = 8;
    arguments.kmer_size = 4;
    arguments.errors = 1;


    std::vector<size_t> const & window_span_begin{0, 4, 5};

    size_t const pattern_begin = 1;
    valik::pattern_bounds expected{};

    // C|GCAAAACG|CGGC
    //
    // 0: [C|GCAAAAC]G |CGGC    don't consider
    // 1: C |[GCAAAACG]|CGGC    consider
    // 2: C |G[CAAAACG |C]GGC   don't consider
    // 3: C |GC[AAAACG |CG]GC   don't consider
    // 4: C |GCA[AAACG |CGG]C   don't consider
    // 5: C |GCAA[AACG |CGGC]   don't consider
    //
    // minimisers[0, 1)
    expected.begin_position = 0;
    expected.end_position = 1;

    valik::threshold threshold_data = valik::make_threshold_data(arguments);
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, threshold_data);

    EXPECT_EQ(bounds.begin_position, expected.begin_position);
    EXPECT_EQ(bounds.end_position, expected.end_position);
}

// ====================================================================================
// seq = CCACGTCGAAGGTT     minimisers = [ACGT, CGAA, AAGG, aacc]
// 	w = 8
// 	k = 4
//
// minimiser	             span_begin     minimiser belongs to windows starting at
// ACGT		CCACGTCGAA           0                           0, 1, 2
// CGAA        CGTCGAAG          3                              3
// AAGG         GTCGAAGGT        4                             4, 5
// aacc           CGAAGGTT       6                              6
// ====================================================================================

TEST(make_pattern_bounds, same_minimiser_consecutive_windows_end)
{
    // special case where pattern ends after the first of multiple consecutive minimisers
    valik::search_arguments arguments{};
    arguments.pattern_size = 12;
    arguments.window_size = 8;
    arguments.kmer_size = 4;
    arguments.errors = 1;

    std::vector<size_t> const & window_span_begin{0, 3, 4, 6};

    size_t const pattern_begin = 0;
    valik::pattern_bounds expected{};

    // |CCACGTCGAAGG|TT
    //
    // 0: |[CCACGTCG]AAGG|TT    consider
    // 1: |C[CACGTCGA]AGG|TT    consider
    // 2: |CC[ACGTCGAA]GG|TT    consider
    // 3: |CCA[CGTCGAAG]G|TT    consider
    // 4: |CCAC[GTCGAAGG]|TT    consider
    // 5: |CCACG[TCGAAGG |T]T   don't consider
    // 6: |CCACGT[CGAAGG |TT]   don't consider
    //
    // minimisers[0, 3)
    expected.begin_position = 0;
    expected.end_position = 3;

    valik::threshold threshold_data = valik::make_threshold_data(arguments);
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, threshold_data);

    EXPECT_EQ(bounds.begin_position, expected.begin_position);
    EXPECT_EQ(bounds.end_position, expected.end_position);
}

// ====================================================================================
// seq = CCACGTCGAAGGTT     minimisers = [ACGT, CGAA, AAGG, aacc]
// 	w = 7
// 	k = 4
//
// minimiser	             span_begin     minimiser belongs to windows starting at
// ACGT		CCACGTCGA            0                          0, 1, 2
// CGAA        CGTCGAAG          3                            3, 4
// AAGG          TCGAAGGT        5                            5, 6
// aacc            GAAGGTT       7                              7
// ====================================================================================

TEST(make_pattern_bounds, odd_lengths)
{
    // normal case where kmer size is not divisible by window or pattern size
    valik::search_arguments arguments{};
    arguments.pattern_size = 11;
    arguments.window_size = 7;
    arguments.kmer_size = 4;
    arguments.errors = 1;

    std::vector<size_t> const & window_span_begin{0, 3, 5, 7};

    size_t const pattern_begin = 2;
    valik::pattern_bounds expected{};

    // CC|ACGTCGAAGGT|T
    //
    // 0: [CC|ACGTC]GAAGGT |T   don't consider
    // 1: C[C|ACGTCG]AAGGT |T   don't consider
    // 2: CC |[ACGTCGA]AGGT|T   consider
    // 3: CC |A[CGTCGAA]GGT|T   consider
    // 4: CC |AC[GTCGAAG]GT|T   consider
    // 5: CC |ACG[TCGAAGG]T|T   consider
    // 6: CC |ACGT[CGAAGGT]|T   consider
    // 7: CC |ACGTC[GAAGGT |T]  don't consider
    //
    // minimisers[0, 3)
    expected.begin_position = 0;
    expected.end_position = 3;

    valik::threshold threshold_data = valik::make_threshold_data(arguments);
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, threshold_data);

    EXPECT_EQ(bounds.begin_position, expected.begin_position);
    EXPECT_EQ(bounds.end_position, expected.end_position);
}
