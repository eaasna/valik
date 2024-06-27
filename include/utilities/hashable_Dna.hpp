#include <cstddef> // for size_t

#include <seqan/basic.h> // for seqan2::Dna

#include <seqan3/alphabet/concept.hpp> // for seqan3::alphabet

template <>
struct seqan3::custom::alphabet<seqan2::Dna>
{
    using alphabet_t = seqan2::Dna;

    constexpr static alphabet_t A_Dna('A');
    constexpr static alphabet_t C_Dna('C');
    constexpr static alphabet_t G_Dna('G');
    constexpr static alphabet_t T_Dna('T');
 
    static constexpr size_t alphabet_size = 4;
 
    static size_t to_rank(alphabet_t const a) noexcept
    {
        return static_cast<size_t>(a);
    }
 
    static constexpr alphabet_t & assign_rank_to(size_t const r, alphabet_t & a) noexcept
    {
        switch (r)
        {
        case 1:
            a = C_Dna;
            return a;
        case 2:
            a = G_Dna;
            return a;
        case 3:
            a = T_Dna;
            return a;
        default:
            a = A_Dna;
            return a;
        }
    }
 
    static constexpr char to_char(alphabet_t const a) noexcept
    {
        switch (a)
        {
        case C_Dna:
            return '1';
        case G_Dna:
            return '2';
        case T_Dna:
            return '3';
        default:
            return '0';
        }
    }
 
    static constexpr alphabet_t & assign_char_to(char const c, alphabet_t & a) noexcept
    {
        switch (c)
        {
        case '1':
            a = C_Dna;
            return a;
        case '2':
            a = G_Dna;
            return a;
        case '3':
            a = T_Dna;
            return a;
        default:
            a = A_Dna;
            return a;
        }
    }
};
 
static_assert(seqan3::alphabet<seqan2::Dna>);
