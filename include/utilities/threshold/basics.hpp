#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <vector>
#include <utility>
#include <filesystem>
#include <map>

#include <cereal/archives/binary.hpp> 
#include <cereal/types/tuple.hpp>

namespace valik
{

constexpr uint8_t query_every = 2; // query every 2nd pattern by default

/**
 *  @brief Assuming k-mers are distributed uniformly randomly, what is the expected count of a k-mer in a bin.
 * 
 * @param bin_size Size of reference segment in bases.
 * @param kmer_size K-mer size.
 * @param information_content How much of the sequence is unique across the genome.
 * 
*/
template <typename var_t, typename par_t>
double expected_kmer_occurrences(var_t const & bin_size,
                                par_t const & kmer_size, 
                                double const information_content = 0.75)
{
    constexpr uint8_t alphabet_size{4};
    return (double) (bin_size - kmer_size + 1) / (double) (pow(alphabet_size, (uint8_t) std::round(kmer_size * information_content)));
}

/**
 * @brief K-mer lemma threshold.
*/
inline size_t kmer_lemma_threshold(size_t const l, uint8_t const k, uint8_t const e)
{
    if (l < k)
        return 0;
    if ((l - k + 1) <= (size_t) e*k)
        return 0;
    
    return l - k + 1 - e * k;
}

/**
 * @brief Definition of the search space for the parameter tuning algorithm.
 * 
 * @param max_errors Maximum number of errors.
 * @param max_thresh Maximum thershold for the number of k-mers that are shared after introducing errors.
 * @param max_len The highest minimum local match length to consider.
 * @param kmer_range The range of k-mer sizes.
*/
struct param_space
{
    constexpr static uint8_t max_errors{15};    
    uint16_t max_thresh{20};
    constexpr static size_t max_len{150};
    constexpr static std::pair<uint8_t, uint8_t> kmer_range{9, 23};

    param_space() noexcept = default;
    param_space(param_space const &) noexcept = default;
    param_space & operator=(param_space const &) noexcept = default;
    param_space & operator=(param_space &&) noexcept = default;
    ~param_space() noexcept = default;

    param_space(uint16_t thresh) : max_thresh(thresh) {}

    uint8_t min_k() const
    {
        return std::get<0>(kmer_range);
    }

    uint8_t max_k() const
    {
        return std::get<1>(kmer_range);
    }

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(max_thresh);
    }
};

/**
 * @brief Number of combinations of n take k.
 * 
 * Same as the number of different error configurations if n=sequence length and k=error_count. 
*/
inline uint64_t combinations(size_t const k, size_t const n)
{
    if (n >= k)
    {
        std::vector<size_t> divisors(k);
        std::iota(divisors.begin(), divisors.end(), 1); // fill vector with 1, 2, ..., k
        std::vector<bool> was_divided(k);

        uint64_t combinations{1};
        for (size_t i = n - k + 1; i <= n; i++)
        {
            // avoid overflowing 64bits
            if (combinations >= (std::numeric_limits<uint64_t>::max() / n))
            {
                for (size_t i{0}; i < divisors.size(); i++)
                {
                    if (!was_divided[i])
                    {
                        if (combinations % divisors[i] == 0)
                        {
                            combinations = combinations / divisors[i];
                            was_divided[i] = true;
                        }
                    }
                }
            }
            combinations *= i;
        }
        
        if (std::find(was_divided.begin(), was_divided.end(), false) != was_divided.end())
        {
            for (size_t i{0}; i < divisors.size(); i++)
            {
                if (!was_divided[i])
                    combinations = combinations / divisors[i];
            }
        }
        return combinations;   // same as (uint64_t) (factorial(n) / (factorial(n - k) * factorial(k)));
    }
    else
        return 0;
}

enum class search_kind {LEMMA, HEURISTIC, MINIMISER, STELLAR};

}   //namespace valik
