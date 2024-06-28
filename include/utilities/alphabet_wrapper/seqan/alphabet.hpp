// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides adaptions for the seqan2 alphabet interface to map seqan3 alphabets.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>
#include <string_view>
#include <vector>

#include <seqan/basic.h>

#include <seqan3/core/platform.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/utility/math.hpp>

#include <utilities/alphabet_wrapper/std/tag_invoke.hpp>
namespace seqan2
{
    template <seqan3::alphabet alphabet_t>
    struct alphabet_adaptor
    {
        alphabet_t _symbol{};
        // make this applicable to simple type adapter of seqan?
        constexpr alphabet_adaptor() = default;

        template <typename char_t>
            requires std::same_as<seqan3::alphabet_char_t<alphabet_t>, char_t> &&
                     requires { seqan3::assign_char_to(char_t{}, alphabet_t{}); }
        constexpr explicit alphabet_adaptor(char_t c) : _symbol{seqan3::assign_char_to(c, alphabet_t{})}
        {}

        template <std::integral rank_t>
            requires (!std::same_as<seqan3::alphabet_char_t<alphabet_t>, rank_t>) &&
                     std::convertible_to<rank_t, seqan3::alphabet_rank_t<alphabet_t>> &&
                     requires { seqan3::assign_rank_to(static_cast<seqan3::alphabet_rank_t<alphabet_t>>(rank_t{}), alphabet_t{}); }
        constexpr explicit alphabet_adaptor(rank_t r) : _symbol{seqan3::assign_rank_to(static_cast<seqan3::alphabet_rank_t<alphabet_t>>(r), alphabet_t{})}
        {}

        template <seqan3::alphabet other_alphabet_t>
            requires (!std::same_as<other_alphabet_t, alphabet_t>) &&
                     seqan3::explicitly_convertible_to<other_alphabet_t, alphabet_t>
        constexpr explicit alphabet_adaptor(alphabet_adaptor<other_alphabet_t> other) noexcept :
            _symbol{static_cast<alphabet_t>(other._symbol)}
        {}

        constexpr bool operator==(alphabet_adaptor const &) const noexcept = default;
        constexpr std::strong_ordering operator<=>(alphabet_adaptor const & other) const noexcept
        {
            return seqan3::to_rank(_symbol) <=> seqan3::to_rank(other._symbol);
        }

        // impicit conversion enough for assignment in seqan2 simple types?
        template <std::integral int_t>
        constexpr operator int_t() const noexcept
        {
            return seqan3::to_rank(_symbol);
        }

        constexpr operator char() const noexcept
        {
            return seqan3::to_char(_symbol);
        }

        constexpr operator alphabet_t() const noexcept
        {
            return _symbol;
        }
    };

    template <seqan3::cereal_output_archive archive_t, typename alphabet_t>
    seqan3::alphabet_rank_t<alphabet_adaptor<alphabet_t>>
    CEREAL_SAVE_MINIMAL_FUNCTION_NAME(archive_t const &, alphabet_adaptor<alphabet_t> const & l)
    {
        return seqan3::to_rank(l);
    }

    template <seqan3::cereal_input_archive archive_t, typename alphabet_t>
    void CEREAL_LOAD_MINIMAL_FUNCTION_NAME(archive_t const &,
                                           alphabet_adaptor<alphabet_t> & l,
                                           seqan3::alphabet_rank_t<alphabet_adaptor<alphabet_t>> const & r)
    {
        seqan3::assign_rank_to(r, l);
    }

    template <typename alphabet_t>
    struct ValueSize<alphabet_adaptor<alphabet_t>>
    {
        using Type = decltype(seqan3::alphabet_size<alphabet_t>);
        static constexpr Type VALUE = seqan3::alphabet_size<alphabet_t>;
    };

    template <typename alphabet_t>
    struct BitsPerValue<alphabet_adaptor<alphabet_t>>
    {
        using Type = decltype(seqan3::alphabet_size<alphabet_t>);
        static constexpr Type VALUE = seqan3::detail::ceil_log2(seqan3::alphabet_size<alphabet_t>);
    };

} //namespace seqan2

namespace jst::contrib
{
    using dna4 =  seqan2::alphabet_adaptor<seqan3::dna4>;  //seqan::Dna;
    using dna5 =  seqan2::alphabet_adaptor<seqan3::dna5>;  //seqan::Dna5;
    using dna15 = seqan2::alphabet_adaptor<seqan3::dna15>; //seqan::Iupac;

    inline dna4 operator""_dna4(char const c) noexcept
    {
        return dna4{c};
    }

    inline std::vector<dna4> operator""_dna4(char const * s, std::size_t n)
    {
        std::vector<dna4> r;
        r.reserve(n);
        std::ranges::copy(std::string_view{s, n} | std::views::transform([] (char c) {
            return dna4{c};
        }), std::back_inserter(r));
        return r;
    }

    inline dna5 operator""_dna5(char const c) noexcept
    {
        return dna5{c};
    }

    inline std::vector<dna5> operator""_dna5(char const * s, std::size_t n)
    {
        std::vector<dna5> r;
        r.reserve(n);
        std::ranges::copy(std::string_view{s, n} | std::views::transform([] (char c) {
            return dna5{c};
        }), std::back_inserter(r));
        return r;
    }
} // namespace jst::contrib

namespace seqan2
{
    template <typename char_t, typename alphabet_value_t, typename spec_t>
    inline seqan3::debug_stream_type<char_t> & operator<<(seqan3::debug_stream_type<char_t> & stream,
                                                          SimpleType<alphabet_value_t, spec_t> const & symbol)
    {
        return stream << seqan3::to_char(symbol);
    }
} // namespace seqan2
//     // use a simple adaptor!
//     // Idea: define a wrapper for seqan3::dna5 as seqan::SimpleType<seqan3::dna5, struct _>
//     // then provide
//     template <seqan3::alphabet alphabet_t>
//     struct adapted;

//     template <seqan3::alphabet alphabet_t>
//     using alphabet_adapter_t = SimpleType<seqan3::alphabet_rank_t<alphabet_t>, adapted<alphabet_t>>;

//     template <typename rank_t, seqan3::alphabet alphabet1_t, seqan3::alphabet alphabet2_t>
//     inline void assign(SimpleType<rank_t, adapted<alphabet1_t>> & target, alphabet2_t & source)
//     {
//         target.value = seqan3::to_rank((alphabet2_t &&) source);
//     }

//     template <typename rank_t, seqan3::alphabet alphabet1_t, seqan3::alphabet alphabet2_t>
//     inline void assign(SimpleType<rank_t, adapted<alphabet1_t>> & target, alphabet2_t const & source)
//     {
//         target.value = seqan3::to_rank((alphabet2_t &&) source);
//     }

//     template <typename rank_t, seqan3::alphabet alphabet_t>
//     inline void assign(char & c_target, SimpleType<rank_t, adapted<alphabet_t>> const & source)
//     {
//         c_target = seqan3::to_char(source);
//     }

//     // need a wrapper that implements seqan::FiniteOrderedAlphabetConcept

//     template <typename rank_t, seqan3::alphabet alphabet_t>
//     struct ValueSize<SimpleType<rank_t, adapted<alphabet_t>>>
//     {
//         using Type = seqan3::detail::min_viable_uint_t<seqan3::alphabet_size<alphabet_t>>;
//         static const Type VALUE = seqan3::alphabet_size<alphabet_t>;
//     };

//     template <typename rank_t, seqan3::alphabet alphabet_t>
//     struct BitsPerValue<SimpleType<rank_t, adapted<alphabet_t>>>
//     {
//         static constexpr size_t bits_per_letter = seqan3::detail::ceil_log2(seqan3::alphabet_size<alphabet_t>);

//         using Type = seqan3::detail::min_viable_uint_t<bits_per_letter>;
//         static const Type VALUE = bits_per_letter;
//     };

//     // sameType(ordValue(val), size);
//     // sameType(valueSize<TValue>(), size);
//     // sameType(ValueSize<TValue>::VALUE, size);

//     // convert integer to alphabet value
//     // val = 0;
//     // val = size;
//     // TValue val2(0);
//     // TValue val3(size);
//     // assign(val, val2);



//     template <seqan3::semialphabet target_t, seqan3::semialphabet source_t>
//         requires std::constructible_from<seqan3::alphabet_rank_t<std::remove_cvref_t<target_t>>, seqan3::alphabet_rank_t<source_t>>
//     struct Convert<target_t, source_t>
//     {
//         using Type = seqan3::alphabet_rank_t<target_t>;
//     };
//     // template <typename TTarget, typename T, typename TSourceValue, typename TSourceSpec>
//     // inline typename RemoveConst_<TTarget>::Type
//     // convertImpl(Convert<TTarget, T> const,
//     //             SimpleType<TSourceValue, TSourceSpec> const & source_)
//     // {
//     //     typename RemoveConst_<TTarget>::Type target_;
//     //     assign(target_, source_);
//     //     return target_;
//     // }
//     // we do not have implicit conversion
//     // so what are we doing now?
//     template <seqan3::semialphabet target_t, typename _t, seqan3::semialphabet source_t>
//         requires (std::constructible_from<seqan3::alphabet_rank_t<std::remove_cvref_t<target_t>>, seqan3::alphabet_rank_t<source_t>>)
//     inline typename Convert<target_t, source_t>::Type
//     convertImpl(Convert<target_t, _t> const, source_t const & source)
//     {
//         return seqan3::to_rank(seqan3::assign_rank_to(seqan3::to_rank(source), target_t{}));
//     }
// }  // namespace seqan

namespace seqan3::custom
{
    template <typename char_t, typename spec_t>
    struct alphabet<seqan2::SimpleType<char_t, spec_t>>
    {
        using alphabet_t = seqan2::SimpleType<char_t, spec_t>;
        using rank_t = typename seqan2::ValueSize<alphabet_t>::Type;

        static constexpr rank_t alphabet_size = seqan2::ValueSize<alphabet_t>::VALUE;

        static rank_t to_rank(alphabet_t const a) noexcept
        {
            return static_cast<rank_t>(a);
        }

        static alphabet_t & assign_rank_to(rank_t const r, alphabet_t & a) noexcept
        {
            a = r; // simple assignment? or convert? or something else?
            return a;
        }

        static char to_char(alphabet_t const a) noexcept
        {
            return static_cast<char>(a);
        }

        static alphabet_t & assign_char_to(char const c, alphabet_t & a) noexcept
        {
            a = c;
            return a;
        }
    };

    template <seqan3::alphabet wrapped_t>
    struct alphabet<seqan2::alphabet_adaptor<wrapped_t>>
    {
        using alphabet_t = seqan2::alphabet_adaptor<wrapped_t>;
        using rank_t = seqan3::alphabet_rank_t<wrapped_t>;
        using char_t = seqan3::alphabet_char_t<wrapped_t>;

        static constexpr rank_t alphabet_size = seqan3::alphabet_size<wrapped_t>;

        constexpr static rank_t to_rank(alphabet_t const a) noexcept
        {
            return seqan3::to_rank(a._symbol);
        }

        constexpr static alphabet_t & assign_rank_to(rank_t const r, alphabet_t & a) noexcept
        {
            seqan3::assign_rank_to(r, a._symbol);
            return a;
        }

        constexpr static char to_char(alphabet_t const a) noexcept
        {
            return seqan3::to_char(a._symbol);
        }

        constexpr static alphabet_t & assign_char_to(char const c, alphabet_t & a) noexcept
        {
            seqan3::assign_char_to(c, a._symbol);
            return a;
        }

        template <seqan3::nucleotide_alphabet _wrapped_t = wrapped_t>
            requires std::same_as<_wrapped_t, wrapped_t>
        constexpr static auto complement(alphabet_t const a)
            noexcept(std::is_nothrow_invocable_v<std::tag_t<seqan3::complement>, wrapped_t const &>)
            -> std::invoke_result_t<std::tag_t<seqan3::complement>, wrapped_t const &>
        {
            return seqan3::complement(a._symbol);
        }

        constexpr static auto char_is_valid(char_t const c) noexcept
        {
            return seqan3::char_is_valid_for<wrapped_t>(c);
        }
    };
} // namespace seqan3::custom

namespace seqan2
{
    template <seqan3::cereal_output_archive archive_t, typename char_t, typename spec_t>
    seqan3::alphabet_rank_t<seqan2::SimpleType<char_t, spec_t>>
    CEREAL_SAVE_MINIMAL_FUNCTION_NAME(archive_t const &, seqan2::SimpleType<char_t, spec_t> const & l)
    {
        return seqan3::to_rank(l);
    }

    template <seqan3::cereal_input_archive archive_t, typename char_t, typename spec_t>
    void CEREAL_LOAD_MINIMAL_FUNCTION_NAME(archive_t const &,
                                           seqan2::SimpleType<char_t, spec_t> & l,
                                           seqan3::alphabet_rank_t<seqan2::SimpleType<char_t, spec_t>> const & r)
    {
        seqan3::assign_rank_to(r, l);
    }
} // namespace seqan2
