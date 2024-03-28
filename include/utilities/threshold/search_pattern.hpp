#pragma once

#include <utilities/threshold/basics.hpp>

#include <cereal/archives/binary.hpp> 

namespace valik
{

/**
 * @brief Define an approximate local alignment of a minimum length and maximum number of errors. 
 *        Same as Stellar epsilon matches.  
*/
struct search_pattern
{
    uint8_t e;
    size_t l;

    search_pattern() noexcept = default;
    search_pattern(search_pattern const &) noexcept = default;
    search_pattern & operator=(search_pattern const &) noexcept = default;
    search_pattern & operator=(search_pattern &&) noexcept = default;
    ~search_pattern() noexcept = default;

    search_pattern(uint8_t const errors, size_t const min_len) : e(errors), l(min_len) {}

    /**
     * @brief Total number of error configurations.
    */
    uint64_t total_conf_count() const
    {
        return combinations(e, l);
    }

    bool operator==(search_pattern const & other) const
    {
        return ((e == other.e) && (l == other.l));
    }

    template <class Archive>
    void serialize(Archive & archive)
    {
        archive(e, l);
    }
};

}   // namespace valik
