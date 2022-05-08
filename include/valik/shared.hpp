#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/threshold/threshold_parameters.hpp>

#include <valik/search/sync_out.hpp>

namespace valik
{

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
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
    std::filesystem::path ref_file{};
    std::filesystem::path ref_out{"reference_metadata.txt"};
    std::filesystem::path seg_out{"reference_segments.txt"};

    size_t overlap{150};
    size_t bins{64};
};

struct build_arguments
{
    uint32_t window_size{23u};
    uint8_t kmer_size{20};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_weight{shape.count()};
    uint8_t threads{1u};

    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    std::filesystem::path out_path{"./"};
    std::string size{};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    bool compressed{false};

    bool from_segments{false};
    std::filesystem::path seg_path{};
    std::filesystem::path ref_meta_path{};
};

struct search_arguments
{
    uint32_t window_size{23u};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_size{shape.size()};
    uint8_t shape_weight{shape.count()};
    uint64_t overlap{};

    uint8_t threads{1u};

    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path query_file{};
    std::filesystem::path index_file{};
    std::filesystem::path out_file{"search.out"};
    std::filesystem::path bin_query_file{};

    // Related to thresholding
    double tau{0.9999};
    double threshold{std::numeric_limits<double>::quiet_NaN()};
    double p_max{0.15};
    double fpr{0.05};
    uint8_t errors{0};
    uint64_t pattern_size{};
    bool treshold_was_set{false};
    bool cache_thresholds{false};

    bool compressed{false};
    bool write_time{false};

    raptor::threshold::threshold_parameters make_threshold_parameters() const noexcept
    {
        return
        {
            .window_size{window_size},
            .shape{shape},
            .pattern_size{pattern_size},
            .errors{errors},
            .percentage{threshold},
            .p_max{p_max},
            .fpr{fpr},
            .tau{tau},
            .cache_thresholds{cache_thresholds},
            .output_directory{index_file.parent_path()}
        };
    }
};

} // namespace valik
