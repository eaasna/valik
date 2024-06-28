// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides container adapter for the seqan container interface.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <optional>
#include <ranges>
#include <type_traits>

#include <seqan/sequence.h>

namespace jst::contrib
{
    template <typename range_t>
    class seqan_container_adapter
    {
    private:
        std::optional<range_t> _range{std::nullopt};

    public:
        using value_type = std::ranges::range_value_t<range_t>;
        using reference = std::ranges::range_reference_t<range_t>;
        using const_reference = std::ranges::range_reference_t<std::remove_const_t<range_t> const>;
        using iterator = std::ranges::iterator_t<range_t>;
        using const_iterator = std::ranges::iterator_t<std::remove_const_t<range_t> const>;
        using size_type = std::make_unsigned_t<std::iter_difference_t<iterator>>;
        using difference_type = std::iter_difference_t<iterator>;

        seqan_container_adapter() = default;
        explicit seqan_container_adapter(range_t range) noexcept(std::is_nothrow_move_constructible_v<range_t>) :
            _range{std::move(range)}
        {
        }

        seqan_container_adapter(seqan_container_adapter const &) = default;
        seqan_container_adapter(seqan_container_adapter &&) = default;
        seqan_container_adapter & operator=(seqan_container_adapter const &) = default;
        seqan_container_adapter & operator=(seqan_container_adapter &&) = default;

        constexpr reference operator[](difference_type const index) noexcept
            requires std::ranges::random_access_range<range_t>
        {
            return (*_range)[index];
        }

        constexpr const_reference operator[](difference_type const index) const noexcept
            requires std::ranges::random_access_range<range_t const>
        {
            return (*_range)[index];
        }

        constexpr iterator begin() noexcept(noexcept(std::ranges::begin(*_range)))
        {
            return std::ranges::begin(*_range);
        }

        constexpr const_iterator begin() const noexcept(noexcept(std::ranges::begin(*_range)))
        {
            return std::ranges::begin(*_range);
        }

        constexpr iterator end() noexcept(noexcept(std::ranges::end(*_range)))
        {
            return std::ranges::end(*_range);
        }

        constexpr const_iterator end() const noexcept(noexcept(std::ranges::end(*_range)))
        {
            return std::ranges::end(*_range);
        }

        constexpr size_type size() const noexcept(noexcept(std::ranges::distance(*_range)))
        {
            return std::ranges::distance(*_range);
        }

        constexpr bool empty() const noexcept(noexcept(std::ranges::empty(*_range)))
        {
            return std::ranges::empty(*_range);
        }
    };

    template <std::ranges::view range_t>
        requires (std::ranges::common_range<range_t> && std::ranges::random_access_range<range_t>)
    constexpr auto make_seqan_container(range_t range) noexcept
    {
        return seqan_container_adapter<std::remove_reference_t<range_t>>{std::move(range)};
    }

    template <std::ranges::view range_t>
    using seqan_container_t = decltype(make_seqan_container(std::declval<range_t>()));

    template <typename range_t>
    inline void
    assign(seqan_container_adapter<range_t> & me, seqan_container_adapter<range_t> const & source)
    {
        seqan2::assign(me, source);
    }

}  // namespace jst::contrib

namespace seqan2
{
    namespace detail
    {
        template <typename range_t>
        using container_t = jst::contrib::seqan_container_adapter<range_t>;
    } // namespace detail

    template <typename range_t>
    struct Value<jst::contrib::seqan_container_adapter<range_t>>
    {
        using Type = typename detail::container_t<range_t>::value_type;
    };
    template <typename range_t>
    struct Value<jst::contrib::seqan_container_adapter<range_t> const>
    {
        using Type = typename detail::container_t<range_t>::value_type;
    };
    template <typename range_t>
    struct Reference<jst::contrib::seqan_container_adapter<range_t>>
    {
        using Type = typename detail::container_t<range_t>::reference;
    };
    template <typename range_t>
    struct Reference<jst::contrib::seqan_container_adapter<range_t> const>
    {
        using Type = typename detail::container_t<range_t>::const_reference;
    };
    template <typename range_t>
    struct GetValue<jst::contrib::seqan_container_adapter<range_t>>
    {
        using Type = typename detail::container_t<range_t>::const_reference;
    };
    template <typename range_t>
    struct GetValue<jst::contrib::seqan_container_adapter<range_t> const>
    {
        using Type = typename detail::container_t<range_t>::const_reference;
    };
    template <typename range_t>
    struct Position<jst::contrib::seqan_container_adapter<range_t>>
    {
        using Type = typename detail::container_t<range_t>::difference_type;
    };
    template <typename range_t>
    struct Position<jst::contrib::seqan_container_adapter<range_t> const>
    {
        using Type = typename detail::container_t<range_t>::difference_type;
    };
    template <typename range_t>
    struct Size<jst::contrib::seqan_container_adapter<range_t>>
    {
        using Type = typename detail::container_t<range_t>::size_type;
    };
    template <typename range_t>
    struct Size<jst::contrib::seqan_container_adapter<range_t> const>
    {
        using Type = typename detail::container_t<range_t>::size_type;
    };
    template <typename range_t>
    struct Iterator<jst::contrib::seqan_container_adapter<range_t>, Standard>
    {
        typedef Iter<detail::container_t<range_t>, StdIteratorAdaptor> Type;
    };
    template <typename range_t>
    struct Iterator<jst::contrib::seqan_container_adapter<range_t> const, Standard>
    {
        typedef Iter<detail::container_t<range_t> const, StdIteratorAdaptor> Type;
    };
    template <typename range_t>
    struct Iterator<jst::contrib::seqan_container_adapter<range_t>, Rooted>
    {

        typedef Iter<detail::container_t<range_t>,
                     AdaptorIterator<
                        Iter<
                            detail::container_t<range_t>,
                            StdIteratorAdaptor
                        >
                    >
                > Type;
    };
    template <typename range_t>
    struct Iterator<jst::contrib::seqan_container_adapter<range_t> const, Rooted>
    {
        typedef Iter<detail::container_t<range_t> const,
                     AdaptorIterator<
                        Iter<
                            detail::container_t<range_t> const,
                            StdIteratorAdaptor
                        >
                    >
                > Type;
    };

    template <typename range_t>
    SEQAN_CONCEPT_IMPL((jst::contrib::seqan_container_adapter<range_t>), (StlContainerConcept));

    template <typename range_t>
    SEQAN_CONCEPT_IMPL((jst::contrib::seqan_container_adapter<range_t> const), (StlContainerConcept));

} // namespace seqan2
