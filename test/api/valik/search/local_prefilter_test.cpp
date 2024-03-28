#include <gtest/gtest.h>

#include <valik/search/local_prefilter.hpp>
#include <valik/shared.hpp>

#include <raptor/threshold/threshold.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

TEST(pattern_begin_positions, read_length_and_pattern_size_are_equal)
{
    // edge case where read_len = pattern_size
    size_t const read_len = 150u;
    uint64_t const pattern_size = read_len;
    uint64_t const query_every = 120u;

    std::vector<size_t> begin_positions{};
    std::vector<size_t> expected{0u}; // pattern interval [0, 150]

    valik::pattern_begin_positions(read_len, pattern_size, query_every, [&](size_t const begin)
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
    uint64_t const query_every = 20u;

    std::vector<size_t> begin_positions{};
    std::vector<size_t> expected = {
        0u, // pattern interval [0, 50]
        0u + query_every /* 20u */, // pattern interval [20, 70]
        20u + query_every /* 40u */, // pattern interval [40, 90]
        60u, // pattern interval [60, 110]
        80u, // pattern interval [80, 130]
        // pattern interval [100, 150]
        100u // ends here, because read_len - pattern_size = 100u
    };

    valik::pattern_begin_positions(read_len, pattern_size, query_every, [&](size_t const begin)
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
    uint64_t const query_every = 30;

    std::vector<size_t> begin_positions{};
    std::vector<size_t> expected = {
        0u, // pattern interval [0, 50]
        0u + query_every /* 30u */, // pattern interval [30, 80]
        30u + query_every /* 60u */, // pattern interval [60, 110]
        90u, // pattern interval [90, 140]
        // pattern interval [100, 150]
        100u // ends here, because read_len - pattern_size = 110u
    };

    valik::pattern_begin_positions(read_len, pattern_size, query_every, [&](size_t const begin)
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
    arguments.shape = seqan3::ungapped{4};
    arguments.errors = 1;

    std::vector<size_t> const & window_span_begin{0, 4, 5};

    size_t const pattern_begin = 0;
    valik::pattern_bounds expected{};


    // |CGCAAAACGCGG|C
    //
    // 0: |[CGC{AAAA}C]GCGG|C     consider (minimiser[0])
    // 1: |C[GC{AAAA}CG]CGG|C     consider
    // 2: |CG[C{AAAA}CGC]GG|C     consider
    // 3: |CGC[{AAAA}CGCG]G|C     consider
    // 4: |CGCA[{AAAC}GCGG]|C     consider (minimiser[1])
    // 5: |CGCAA[{AACG}CGG |C]    don't consider (minimiser[2])
    //
    // minimisers[0, 2)
    expected.begin_position = 0;
    expected.end_position = 2;


    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, thresholder);

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
    arguments.shape = seqan3::ungapped{4};
    arguments.errors = 1;


    std::vector<size_t> const & window_span_begin{0, 4, 5};

    size_t const pattern_begin = 1;
    valik::pattern_bounds expected{};

    // C|GCAAAACGCGGC|
    //
    // 0: [C|GC{AAAA}C]GCGGC |     don't consider (minimiser[0])
    // 1:  C|[GC{AAAA}CG]CGGC|     consider (still minimiser[0])
    // 2:  C|G[C{AAAA}CGC]GGC|     consider
    // 3:  C|GC[{AAAA}CGCG]GC|     consider
    // 4:  C|GCA[{AAAC}GCGG]C|     consider (minimiser[1])
    // 5:  C|GCAA[{AACG}CGGC]|     consider (minimiser[2])
    //
    // minimisers[0, 3)
    expected.begin_position = 0;
    expected.end_position = 3;

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, thresholder);

    EXPECT_EQ(bounds.begin_position, expected.begin_position);
    EXPECT_EQ(bounds.end_position, expected.end_position);
}

TEST(make_pattern_bounds, pattern_equals_window)
{
    // special case where pattern_size == window_size
    valik::search_arguments arguments{};
    arguments.pattern_size = 8;
    arguments.window_size = 8;
    arguments.shape = seqan3::ungapped{4};
    arguments.errors = 1;


    std::vector<size_t> const & window_span_begin{0, 4, 5};

    size_t const pattern_begin = 1;
    valik::pattern_bounds expected{};

    // C|GCAAAACG|CGGC
    //
    // 0: [C|GC{AAAA}C]G |CGGC    don't consider (minimiser[0])
    // 1: C |[GC{AAAA}CG]|CGGC    consider (still minimiser[0])
    // 2: C |G[C{AAAA}CG |C]GGC   don't consider
    // 3: C |GC[{AAAA}CG |CG]GC   don't consider
    // 4: C |GCA[{AAAC}G |CGG]C   don't consider (minimiser[1])
    // 5: C |GCAA[{AACG} |CGGC]   don't consider (minimiser[2])
    //
    // minimisers[0, 1)
    expected.begin_position = 0;
    expected.end_position = 1;

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, thresholder);

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
    arguments.shape = seqan3::ungapped{4};
    arguments.errors = 1;

    std::vector<size_t> const & window_span_begin{0, 3, 4, 6};

    size_t const pattern_begin = 0;
    valik::pattern_bounds expected{};

    // |CCACGTCGAAGG|TT
    //
    // 0: |[CC{ACGT}CG]AAGG|TT    consider (minimiser[0])
    // 1: |C[C{ACGT}CGA]AGG|TT    consider
    // 2: |CC[{ACGT}CGAA]GG|TT    consider
    // 3: |CCA[CGT{CGAA}G]G|TT    consider (minimiser[1])
    // 4: |CCAC[GTCG{AAGG}]|TT    consider (minimiser[2])
    // 5: |CCACG[TCG{AAGG} |T]T   don't consider
    // 6: |CCACGT[CGAA{GG  |TT}]  don't consider (minimiser[3] = aacc)
    //
    // minimisers[0, 3)
    expected.begin_position = 0;
    expected.end_position = 3;

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, thresholder);

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
    arguments.shape = seqan3::ungapped{4};
    arguments.errors = 1;

    std::vector<size_t> const & window_span_begin{0, 3, 5, 7};

    size_t const pattern_begin = 2;
    valik::pattern_bounds expected{};

    // CC|ACGTCGAAGGT|T
    //
    // 0: [CC|{ACGT}C]GAAGGT |T   don't consider (minimiser[0])
    // 1: C[C|{ACGT}CG]AAGGT |T   don't consider
    // 2: CC |[{ACGT}CGA]AGGT|T   consider (still minimiser[0])
    // 3: CC |A[CGT{CGAA}]GGT|T   consider (minimiser[1])
    // 4: CC |AC[GT{CGAA}G]GT|T   consider
    // 5: CC |ACG[TCG{AAGG}]T|T   consider (minimiser[2])
    // 6: CC |ACGT[CG{AAGG}T]|T   consider
    // 7: CC |ACGTC[GAA{GGT  |T}] don't consider (minimiser[3] = aacc)
    //
    // minimisers[0, 3)
    expected.begin_position = 0;
    expected.end_position = 3;

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};
    valik::pattern_bounds bounds = valik::make_pattern_bounds(pattern_begin, arguments, window_span_begin, thresholder);

    EXPECT_EQ(bounds.begin_position, expected.begin_position);
    EXPECT_EQ(bounds.end_position, expected.end_position);
}
