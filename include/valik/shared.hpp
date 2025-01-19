#pragma once

#include <filesystem>
#include <vector>

#include <utilities/threshold/basics.hpp>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/threshold/threshold.hpp>
#include <raptor/threshold/threshold_parameters.hpp>

#include <dream_stellar/shared.hpp>

namespace valik
{

constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

/**
 * @brief Function that rounds the chosen segment count to the closest multiple of 64.
 *
 * @param n Segment count.
 */
constexpr static size_t adjust_bin_count(size_t const & n)
{
    if (n == std::numeric_limits<uint32_t>::max())
        return 64u;
    int remainder = n % 64;
    if (remainder == 0)
        return n;
    else if (remainder < 32)
        return std::max((uint32_t) n - remainder, 128u) - 64;  // previous multiple of 64
    else
        return n + 64 - remainder;  // next multiple of 64
}

//!\brief Strong type for passing the window size.
struct window { uint64_t v; };
//!\brief Strong type for passing the kmer size.
struct kmer { uint8_t v; };
//!\brief Strong type for passing number of bins.
struct bins { uint64_t v; };
//!\brief Strong type for passing number of bits.
struct bits { uint64_t v; };
//!\brief Strong type for passing number of hash functions.
struct hashes { uint64_t v; };

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct split_arguments
{
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path db_file{};
    std::filesystem::path meta_out{"metadata.bin"};

    size_t pattern_size{150};
    uint32_t seg_count{64};
    uint32_t seg_count_in{std::numeric_limits<uint32_t>::max()};
    float error_rate{0.05};
    uint8_t errors{0};
    float fpr{0.05};
    uint8_t kmer_size{std::numeric_limits<uint8_t>::max()};
    std::string shape_str{};
    seqan3::shape shape{};
    uint8_t shape_weight{};
    
    bool metagenome{false};
    std::filesystem::path ref_meta_path{};
    bool write_out{false};
    bool split_query{false};
    bool only_split{false};
    bool verbose{false};
};

struct build_arguments
{
    uint32_t window_size{23u};
    uint8_t kmer_size{20};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_weight{shape.count()};
    uint8_t threads{1u};

    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path out_path;
    std::filesystem::path out_dir{"./"};
    float fpr{0.05};
    std::string size{};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    bool fast{false};
    bool manual_parameters{false};
    bool input_is_minimiser{false};

    uint8_t kmer_count_min_cutoff{2};
    uint8_t kmer_count_max_cutoff{64};
    bool use_filesize_dependent_cutoff{false};

    std::filesystem::path ref_meta_path{};
    bool verbose{false};
};

struct minimiser_threshold_arguments
{
    virtual ~minimiser_threshold_arguments() = 0;   // make an abstract base struct

    double tau{0.9999};
    double p_max{0.15};
    double fpr{0.05};
    uint8_t errors{0};
    size_t pattern_size{100};
    double threshold_percentage{std::numeric_limits<double>::quiet_NaN()};
    bool threshold_was_set{false};
    bool cache_thresholds{false};

    protected:
        // prevent creating, assigning or moving base struct instances
        minimiser_threshold_arguments() = default;
        minimiser_threshold_arguments(minimiser_threshold_arguments const&) = default;
        minimiser_threshold_arguments(minimiser_threshold_arguments&&) = default;
        minimiser_threshold_arguments& operator=(minimiser_threshold_arguments const&) = default;
        minimiser_threshold_arguments& operator=(minimiser_threshold_arguments&&) = default;
};

inline minimiser_threshold_arguments::~minimiser_threshold_arguments() = default;

struct search_profile_arguments
{
    virtual ~search_profile_arguments() = 0;   // make an abstract base struct

    //!TODO: deduce this automatically
    bool split_query{false};
    search_kind search_type{search_kind::LEMMA};
    double fnr;

    protected:
        // prevent creating, assigning or moving base struct instances
        search_profile_arguments() = default;
        search_profile_arguments(search_profile_arguments const&) = default;
        search_profile_arguments(search_profile_arguments&&) = default;
        search_profile_arguments& operator=(search_profile_arguments const&) = default;
        search_profile_arguments& operator=(search_profile_arguments&&) = default;
};

inline search_profile_arguments::~search_profile_arguments() = default;

struct search_arguments final : public minimiser_threshold_arguments, search_profile_arguments, public dream_stellar::StellarOptions
{
    ~search_arguments() override = default;
    search_arguments() = default;
    search_arguments(search_arguments const&) = default;
    search_arguments(search_arguments&&) = default;
    search_arguments& operator=(search_arguments const&) = default;
    search_arguments& operator=(search_arguments&&) = default;

    uint32_t window_size{23u};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_size{shape.size()};
    uint8_t shape_weight{shape.count()};
    uint8_t query_every{2};
    size_t threshold{};
    uint32_t seg_count{64};
    uint32_t seg_count_in{std::numeric_limits<uint32_t>::max()};
    uint64_t max_segment_len{(uint64_t) 1e5};

    uint8_t threads{1u};

    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path query_file{};
    std::filesystem::path index_file{};
    std::filesystem::path all_matches{};
    std::filesystem::path out_file{"search.gff"};

    bool write_time{false};
    bool fast{false};
    bool verbose{false};
    bool keep_best_repeats{false};
    double best_bin_entropy_cutoff{0.25};
    bool keep_all_repeats{false};
    bool stellar_only{false};

    size_t cart_max_capacity{1000};
    size_t max_queued_carts{std::numeric_limits<size_t>::max()};

    raptor::threshold::threshold_parameters make_threshold_parameters() const noexcept
    {
        return
        {
            .window_size{window_size},
            .shape{shape},
            .query_length{pattern_size},
            .errors{errors},
            .percentage{threshold_percentage},
            .p_max{p_max},
            .fpr{fpr},
            .tau{tau},
            .cache_thresholds{cache_thresholds},
            .output_directory{index_file.parent_path()}
        };
    }

    float error_rate{};
    std::filesystem::path ref_meta_path{};
    bool distribute{false};
};

} // namespace valik
