// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts and CPOs for the pattern adaption.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utilities/alphabet_wrapper/std/tag_invoke.hpp>

namespace jst::contrib
{

    namespace _set_up
    {
        inline constexpr struct _cpo
        {
            template <typename operation_t, typename finder_t>
                requires std::tag_invocable<_cpo, operation_t &, finder_t &>
            constexpr auto operator()(operation_t & operation, finder_t & finder) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, operation_t &, finder_t &>)
                -> std::tag_invoke_result_t<_cpo, operation_t &, finder_t &>
            {
                return std::tag_invoke(_cpo{}, operation, finder);
            }
        } set_up{};
    } // namespace _set_up
    using _set_up::set_up;

    namespace _tear_down
    {
        inline constexpr struct _cpo
        {
            template <typename operation_t, typename finder_t>
                requires std::tag_invocable<_cpo, operation_t &, finder_t &>
            constexpr auto operator()(operation_t & operation, finder_t & finder) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, operation_t &, finder_t &>)
                -> std::tag_invoke_result_t<_cpo, operation_t &, finder_t &>
            {
                return std::tag_invoke(_cpo{}, operation, finder);
            }
        } tear_down{};
    } // namespace _tear_down
    using _tear_down::tear_down;

}  // namespace jst::contrib
