#pragma once

#include <utilities/threshold/param_set.hpp>
#include <utilities/threshold/search_pattern.hpp>

#include <unordered_map>
#include <cereal/archives/binary.hpp> 
#include <cereal/types/unordered_map.hpp>


namespace valik
{

struct error_threshold
{
    param_set params;
    search_pattern pattern;
    bool is_heuristic{false};
    double fnr;
    double fp_per_pattern;
    uint64_t max_segment_len;

    error_threshold() noexcept = default;
    error_threshold(error_threshold const &) noexcept = default;
    error_threshold & operator=(error_threshold const &) noexcept = default;
    error_threshold & operator=(error_threshold &&) noexcept = default;
    ~error_threshold() noexcept = default;

    error_threshold(param_set const & best_params, search_pattern const & eps_pattern, bool const is_approximate, 
                    double const false_neg, double const false_pos, double const max_len) :
                    params(best_params), pattern(eps_pattern), is_heuristic(is_approximate), fnr(false_neg),
                    fp_per_pattern(false_pos), max_segment_len(max_len) {}

    template <class Archive>
    void serialize(Archive & archive)
    {
        archive(params, pattern, is_heuristic, fnr, fp_per_pattern, max_segment_len);
    }

    void print()
    {
        if (is_heuristic)
            std::cout << "heuristic";
        else
            std::cout << "kmer lemma";
        
        std::cout << '\t' << std::to_string(params.t) << '\t' << fnr << '\t'
                  << fp_per_pattern << '\t' << max_segment_len << '\n';
    }
};

struct kmer_thresholds
{
    size_t k;
    size_t l;
    std::unordered_map<uint8_t, error_threshold> error_table;
    uint8_t max_errors{0};

    kmer_thresholds() noexcept = default;
    kmer_thresholds(kmer_thresholds const &) noexcept = default;
    kmer_thresholds & operator=(kmer_thresholds const &) noexcept = default;
    kmer_thresholds & operator=(kmer_thresholds &&) noexcept = default;
    ~kmer_thresholds() noexcept = default;

    kmer_thresholds(size_t const kmer_size, size_t const min_len) : k(kmer_size), l(min_len) {}

    kmer_thresholds(std::filesystem::path const & filepath)
    {
        load(filepath);
    }

    void add_error_rate(uint8_t const errors, error_threshold const & error_thresh)
    {        
        if (errors > max_errors)
            max_errors = errors;
        if (error_thresh.pattern.l != l)
            throw std::runtime_error("Incompatible pattern length");
        if (error_thresh.params.k != k)
            throw std::runtime_error("Incompatible kmer size");
        error_table.insert({errors, error_thresh});
    }

    error_threshold get_error_thresh(uint8_t const errors)
    {
        if (error_table.find(errors) == error_table.end())
            throw std::runtime_error("Error count " + std::to_string(errors) + " out of precalculated range");

        return error_table.at(errors);
    }

    /**
     * @brief Serialize the kmer thresholds struct.
     *
     * @param filepath Output file path.
     */
    void save(std::filesystem::path const & filepath) const
    {
        std::ofstream os(filepath, std::ios::binary);
        cereal::BinaryOutputArchive archive(os);
        archive(k, l, error_table);
    }
      
    /**
     * @brief Deserialise the kmer thresholds struct.
     *
     * @param filepath Input file path.
     */
    void load(std::filesystem::path const & filepath)
    {
        std::ifstream is(filepath, std::ios::binary);
        cereal::BinaryInputArchive archive(is);
        archive(k, l, error_table);
        max_errors = error_table.size();
    }

    void print()
    {
        std::cout.precision(3);
        std::cout << "\nRecommended shared " << std::to_string(k) << "-mer thresholds for different error rates\n";
        std::cout << "error_rate\tthreshold_kind\tthreshold\tFNR\tFP_per_pattern\tmax_segment_len\n";

        for (uint8_t er{1}; er <= max_errors; er++)
        {
            std::cout << er / (double) l << '\t';
            auto error_thresh = error_table.find(er);
            if (error_thresh != error_table.end())
                error_thresh->second.print();
        }
    }
};

} // namespace valik
