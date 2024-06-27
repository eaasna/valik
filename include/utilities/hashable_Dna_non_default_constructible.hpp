#pragma once

#include <cstddef>     // for size_t
#include <type_traits> // for std::type_identity

#include <seqan/basic.h> // for seqan::Dna
 
#include <seqan3/alphabet/concept.hpp> // for seqan3::alphabet
 
namespace valik
{
 
class hashable_Dna
{
    using alphabet_t = seqan2::Dna;
public:
    uint8_t rank;
 
    hashable_Dna() = delete;
    constexpr hashable_Dna(hashable_Dna const &) = default;
    constexpr hashable_Dna & operator=(hashable_Dna const &) = default;
 
    constexpr hashable_Dna(bool rank) : rank{rank}
    {}
 
    constexpr friend bool operator==(hashable_Dna lhs, hashable_Dna rhs)
    {
        return lhs.rank == rhs.rank;
    }
    constexpr friend bool operator!=(hashable_Dna lhs, hashable_Dna rhs)
    {
        return lhs.rank != rhs.rank;
    }
    constexpr friend bool operator<=(hashable_Dna lhs, hashable_Dna rhs)
    {
        return lhs.rank <= rhs.rank;
    }
    constexpr friend bool operator>=(hashable_Dna lhs, hashable_Dna rhs)
    {
        return lhs.rank >= rhs.rank;
    }
    constexpr friend bool operator<(hashable_Dna lhs, hashable_Dna rhs)
    {
        return lhs.rank < rhs.rank;
    }
    constexpr friend bool operator>(hashable_Dna lhs, hashable_Dna rhs)
    {
        return lhs.rank > rhs.rank;
    }
};
 
constexpr size_t alphabet_size(std::type_identity<hashable_Dna> const &) noexcept
{
    return 4;
}
 
constexpr uint8_t to_rank(hashable_Dna const a) noexcept
{
    return a.rank;
}
 
constexpr hashable_Dna & assign_rank_to(uint8_t const r, hashable_Dna & a) noexcept
{
    a.rank = r;
    return a;
}
 
constexpr char to_char(hashable_Dna const a) noexcept
{
    switch (a.rank)
    {
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            return 'A';
    }
}
 
constexpr hashable_Dna & assign_char_to(char const c, hashable_Dna & a) noexcept
{
    switch (c)
    {
    case 'C':
        a.rank = 1;
        return a;
    case 'G':
        a.rank = 2;
        return a;
    case 'T':
        a.rank = 3;
        return a;
    default:
        a.rank = 0;
        return a;
    }
}
 
constexpr bool char_is_valid_for(char const c, std::type_identity<hashable_Dna> const &) noexcept
{
    switch (c)
    {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
        return true;
    default:
        return false;
    }
}
 
} // namespace valik
 
static_assert(seqan3::alphabet_size<valik::hashable_Dna> == 4);
static_assert(seqan3::char_is_valid_for<valik::hashable_Dna>('T'));
static_assert(!seqan3::char_is_valid_for<valik::hashable_Dna>('!'));
static_assert(seqan3::semialphabet<valik::hashable_Dna>);
static_assert(seqan3::alphabet<valik::hashable_Dna>);
