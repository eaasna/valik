// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <type_traits>
#include <utility>

namespace std
{
    // ----------------------------------------------------------------------------
    // Implementation of tag_invoke
    // ----------------------------------------------------------------------------

    namespace _tag_invoke
    {
        // poison pill
        void tag_invoke();

        struct _fn
        {
            template <typename cpo_t, typename... args_t>
            constexpr auto operator()(cpo_t cpo, args_t &&...args) const
                noexcept(noexcept(tag_invoke(std::forward<cpo_t>(cpo), std::forward<args_t>(args)...)))
                    -> decltype(tag_invoke(std::forward<cpo_t>(cpo), std::forward<args_t>(args)...))
            {
                return tag_invoke((cpo_t &&)cpo, (args_t &&)args...);
            }
        };

        template <typename cpo_t, typename... args_t>
        using tag_invoke_result_t = decltype(tag_invoke(std::declval<cpo_t &&>(), std::declval<args_t &&>()...));

        using yes_type = char;
        using no_type = char (&)[2];

        template <typename cpo_t, typename... args_t>
        auto try_tag_invoke(int) //
            noexcept(noexcept(tag_invoke(std::declval<cpo_t &&>(), std::declval<args_t &&>()...)))
                -> decltype(static_cast<void>(tag_invoke(std::declval<cpo_t &&>(), std::declval<args_t &&>()...)), yes_type{});

        template <typename cpo_t, typename... args_t>
        no_type try_tag_invoke(...) noexcept(false);

        template <template <typename...> class t, typename... args_t>
        struct defer
        {
            using type = t<args_t...>;
        };

        struct empty
        {
        };
    } // namespace _tag_invoke

    // ----------------------------------------------------------------------------
    // Definition of tag_invoke
    // ----------------------------------------------------------------------------

    namespace _tag_invoke_cpo
    {
        inline constexpr _tag_invoke::_fn tag_invoke{};
    }
    using namespace _tag_invoke_cpo;

    // ----------------------------------------------------------------------------
    // Utility functions for tag_invoke
    // ----------------------------------------------------------------------------

    template <auto &cpo_t>
    using tag_t = std::remove_cvref_t<decltype(cpo_t)>;

    using _tag_invoke::tag_invoke_result_t;

    template <typename cpo_t, typename... args_t>
    inline constexpr bool is_tag_invocable_v =
        (sizeof(_tag_invoke::try_tag_invoke<cpo_t, args_t...>(0)) ==
         sizeof(_tag_invoke::yes_type));

    template <typename cpo_t, typename... args_t>
    struct tag_invoke_result
        : std::conditional_t<
              is_tag_invocable_v<cpo_t, args_t...>,
              _tag_invoke::defer<tag_invoke_result_t, cpo_t, args_t...>,
              _tag_invoke::empty>
    {
    };

    template <typename cpo_t, typename... args_t>
    using is_tag_invocable = std::bool_constant<is_tag_invocable_v<cpo_t, args_t...>>;

    template <typename cpo_t, typename... args_t>
    inline constexpr bool is_nothrow_tag_invocable_v =
        noexcept(_tag_invoke::try_tag_invoke<cpo_t, args_t...>(0));

    template <typename cpo_t, typename... args_t>
    using is_nothrow_tag_invocable =
        std::bool_constant<is_nothrow_tag_invocable_v<cpo_t, args_t...>>;

    template <typename cpo_t, typename... args_t>
    concept tag_invocable =
        (sizeof(_tag_invoke::try_tag_invoke<cpo_t, args_t...>(0)) ==
         sizeof(_tag_invoke::yes_type));
} // namespace std
