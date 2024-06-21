#pragma once

#include <utilities/threshold/param_set.hpp>
#include <utilities/threshold/search_pattern.hpp>

#include <unordered_map>
#include <cereal/archives/binary.hpp> 
#include <cereal/types/unordered_map.hpp>

namespace valik
{

/**
 * @brief Search accuracy for an epsilon match definition and chosen parameters k and t.
 * 
 * \param params Chosen kmer size and threshold
 * \param pattern Definition of an epsilon match with a minimum length and number of errors
 * \param search_type One of: LEMMA, HEURISTIC, STELLAR
 * \param fnr False negative rate of search, shared k-mer count is less than threshold although there exists an epsilon match
 * \param fp_per_pattern Probability of a query pattern appearing spuriously in one of the reference bins
 * \param max_segment_len Maximum size of a query segment that does not appear spuriously in more than 10% of reference bins
 * 
*/
struct search_error_profile
{
    param_set params;
    search_pattern pattern;
    search_kind search_type{search_kind::LEMMA};
    double fnr;
    double fp_per_pattern;

    search_error_profile() noexcept = default;
    search_error_profile(search_error_profile const &) noexcept = default;
    search_error_profile & operator=(search_error_profile const &) noexcept = default;
    search_error_profile & operator=(search_error_profile &&) noexcept = default;
    ~search_error_profile() noexcept = default;

    search_error_profile(param_set const & best_params, search_pattern const & eps_pattern, search_kind const which_search, 
                    double const false_neg, double const false_pos) :
                    params(best_params), pattern(eps_pattern), search_type(which_search), fnr(false_neg),
                    fp_per_pattern(false_pos) {}
    
    search_error_profile(param_set const & best_params, search_pattern const & eps_pattern, search_kind const which_search) :
                    params(best_params), pattern(eps_pattern), search_type(which_search) 
    {
        if (which_search != search_kind::STELLAR)
            throw std::runtime_error("Provide FNR, FPR and maximum segment length.");
        fnr = 0.0;
    }

    template <class Archive>
    void serialize(Archive & archive)
    {
        auto type_int = static_cast<int>(search_type);
        archive(params, pattern, type_int, fnr, fp_per_pattern);
        search_type = static_cast<search_kind>(type_int);
    }

    void print() const
    {
        switch (search_type)
        {
            case search_kind::STELLAR: std::cout << "stellar-only" << "\tN/A\t" << fnr << "\tN/A\tN/A\n"; break;
            default:
            {
                switch (search_type)
                {
                    case search_kind::HEURISTIC: std::cout << "heuristic"; break;
                    case search_kind::LEMMA: std::cout << "kmer lemma"; break;
                    default: break;
                }
                std::cout << '\t' << std::to_string(params.t) << '\t' << fnr << '\t' << fp_per_pattern << '\n';
            }
        }
    }
};

/**
 * @brief For a k-mer size and min length gather the search accuracy (FNR, FPR) depending on error rate. 
 * 
 * \param k kmer size
 * \param l minimum local match lenght
 * \param error_table search metrics per error rate
 * \param max_errors maximum error count 
*/
struct search_kmer_profile
{
    size_t k;
    size_t l;
    std::unordered_map<uint8_t, search_error_profile> error_table;
    uint8_t max_errors{0};

    search_kmer_profile() noexcept = default;
    search_kmer_profile(search_kmer_profile const &) noexcept = default;
    search_kmer_profile & operator=(search_kmer_profile const &) noexcept = default;
    search_kmer_profile & operator=(search_kmer_profile &&) noexcept = default;
    ~search_kmer_profile() noexcept = default;

    search_kmer_profile(size_t const kmer_size, size_t const min_len) : k(kmer_size), l(min_len) {}

    search_kmer_profile(std::filesystem::path const & filepath)
    {
        load(filepath);
    }

    void add_error_rate(uint8_t const errors, search_error_profile const & error_thresh)
    {        
        if (errors > max_errors)
            max_errors = errors;
        if (error_thresh.pattern.l != l)
            throw std::runtime_error("Incompatible pattern length");
        if (error_thresh.params.k != k)
            throw std::runtime_error("Incompatible kmer size");
        error_table.insert({errors, error_thresh});
    }

    search_error_profile get_error_profile(uint8_t const errors)
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

    double max_error_rate() const
    {
        return max_errors / (double) l;
    }

    void print() const 
    {
        std::cout.precision(3);
        std::cout << "\nRecommended shared " << std::to_string(k) << "-mer thresholds for matches of (min_length=" << std::to_string(l) 
                  << "; max_error_rate=" << max_error_rate() << ")\n";
        std::cout << "errors\tthreshold_kind\tthreshold\tFNR\tFP_per_pattern\n";

        for (uint8_t er{0}; er <= max_errors; er++)
        {
            std::cout << std::to_string(er) << '\t';
            auto error_thresh = error_table.find(er);
            if (error_thresh != error_table.end())
                error_thresh->second.print();
        }
    }
};

} // namespace valik
