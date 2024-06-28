// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides matcher concepts.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <utilities/alphabet_wrapper/std/tag_invoke.hpp>

namespace jst::contrib
{
    // ----------------------------------------------------------------------------
    // Operation CPOs for matcher
    // ----------------------------------------------------------------------------

    // window_size
    namespace _window_size {
        inline constexpr struct _cpo  {
            template <typename matcher_t>
                requires std::tag_invocable<_cpo, matcher_t>
            constexpr auto operator()(matcher_t && matcher) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, matcher_t>)
                -> std::tag_invoke_result_t<_cpo, matcher_t>
            {
                return std::tag_invoke(_cpo{}, (matcher_t &&)matcher);
            }

        private:

            template <typename matcher_t>
                requires requires (matcher_t && matcher) { { ((matcher_t &&)matcher).window_size() } -> std::integral; }
            constexpr friend auto tag_invoke(_cpo, matcher_t && matcher)
                noexcept(noexcept(std::declval<matcher_t &&>().window_size()))
                -> decltype(std::declval<matcher_t &&>().window_size())
            {
                return ((matcher_t &&)matcher).window_size();
            }
        } window_size;
    } // namespace _window_size
    using _window_size::window_size;

    template <typename matcher_t>
    using window_size_t = std::invoke_result_t<_window_size::_cpo, matcher_t>;

    // capture
    namespace _capture {
        inline constexpr struct _cpo  {
            template <typename matcher_t>
                requires std::tag_invocable<_cpo, matcher_t const &>
            constexpr auto operator()(matcher_t const & matcher) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, matcher_t const &>)
                -> std::tag_invoke_result_t<_cpo, matcher_t const &>
            {
                return std::tag_invoke(_cpo{}, matcher);
            }
        private:

            template <typename matcher_t>
                requires requires (matcher_t && matcher) { { ((matcher_t &&)matcher).capture() }; }
            constexpr friend auto tag_invoke(_cpo, matcher_t && matcher)
                noexcept(noexcept(std::declval<matcher_t &&>().capture()))
                -> decltype(std::declval<matcher_t &&>().capture())
            {
                return ((matcher_t &&)matcher).capture();
            }
        } capture;
    } // namespace _capture
    using _capture::capture;

    // restore
    namespace _restore {
        inline constexpr struct _cpo  {
            template <typename matcher_t, typename state_t>
                requires std::tag_invocable<_cpo, matcher_t &, state_t>
            constexpr auto operator()(matcher_t & matcher, state_t && state) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, matcher_t &, state_t>)
                -> std::tag_invoke_result_t<_cpo, matcher_t &, state_t>
            {
                return std::tag_invoke(_cpo{}, matcher, (state_t &&) state);
            }
        private:

            template <typename matcher_t, typename state_t>
                requires requires (matcher_t && m, state_t && s) { { ((matcher_t &&)m).restore((state_t &&) s) }; }
            constexpr friend void tag_invoke(_cpo, matcher_t && matcher, state_t && state)
                noexcept(noexcept(std::declval<matcher_t &&>().restore(std::declval<state_t &&>())))
            {
                ((matcher_t &&)matcher).restore((state_t &&)state);
            }
        } restore;
    } // namespace _restore
    using _restore::restore;

    template <typename matcher_t>
    using matcher_state_t = std::remove_cvref_t<std::invoke_result_t<_capture::_cpo, matcher_t>>;

    // aggregate
    namespace _aggregate {
        inline constexpr struct _cpo  {
            template <typename state1_t, typename state2_t>
                requires std::tag_invocable<_cpo, state1_t, state2_t>
            constexpr auto operator()(state1_t && state1, state2_t && state2) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, state1_t, state2_t>)
                -> std::tag_invoke_result_t<_cpo, state1_t, state2_t>
            {
                return std::tag_invoke(_cpo{}, (state1_t &&) state1, (state2_t &&) state2);
            }
        } aggregate;
    } // namespace _aggregate
    using _aggregate::aggregate;

    // ----------------------------------------------------------------------------
    // Concept defintions for matcher
    // ----------------------------------------------------------------------------

    template <typename matcher_t>
    concept window_matcher = std::copyable<std::remove_cvref_t<matcher_t>> && requires (matcher_t && matcher)
    {
        { jst::contrib::window_size((matcher_t &&) matcher)} -> std::integral;
    };

    namespace detail {
        template <typename matcher_t>
        concept stateful_matcher = requires
        {
            typename jst::contrib::matcher_state_t<matcher_t>;
            requires std::semiregular<jst::contrib::matcher_state_t<matcher_t>>;
        };
    } // namespace detail

    template <typename matcher_t>
    concept restorable_matcher = window_matcher<matcher_t> &&
                                 detail::stateful_matcher<matcher_t> &&
                                 requires (matcher_t && matcher)
    {
        { jst::contrib::capture((matcher_t &&) matcher) } -> std::convertible_to<jst::contrib::matcher_state_t<matcher_t>>;
        { jst::contrib::restore(matcher, std::declval<jst::contrib::matcher_state_t<matcher_t>>()) };
    };

    template <typename t1, typename t2>
    concept reducable_with = std::common_with<std::remove_cvref_t<t1>, std::remove_cvref_t<t2>> &&
                             requires (std::remove_cvref_t<t1> const & lhs, std::remove_cvref_t<t2> const & rhs)
    {
        { jst::contrib::aggregate(lhs, rhs) } -> std::convertible_to<std::common_type_t<std::remove_cvref_t<t1>, std::remove_cvref_t<t2>>>;
    };

    template <typename state_t>
    concept reducable_state = reducable_with<state_t, state_t>;

    template <typename matcher_t, typename ...args_t>
    concept online_matcher_for = window_matcher<matcher_t> && std::invocable<matcher_t, args_t...>;
}  // namespace jst::contrib
