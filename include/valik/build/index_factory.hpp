#pragma once

#include <valik/build/call_parallel_on_bins.hpp>
#include <valik/index.hpp>
#include <valik/split/metadata.hpp>

#include <seqan3/search/views/minimiser_hash.hpp>

namespace valik
{

template <bool compressed>
class index_factory
{
public:
    index_factory() = default;
    index_factory(index_factory const &) = default;
    index_factory(index_factory &&) = default;
    index_factory & operator=(index_factory const &) = default;
    index_factory & operator=(index_factory &&) = default;
    ~index_factory() = default;

    explicit index_factory(build_arguments const & args) : arguments{std::addressof(args)} {}

    template <typename view_t = int>
    [[nodiscard]] auto operator()(view_t && hash_filter_view = 0) const
    {
        auto tmp = construct(std::move(hash_filter_view));

        if constexpr (!compressed)
            return tmp;
        else
            return valik_index<index_structure::ibf_compressed>{std::move(tmp)};
    }

private:
    build_arguments const * const arguments{nullptr};

    template <typename view_t>
    auto construct(view_t && hash_filter_view) const
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        assert(arguments != nullptr);

        valik_index<> index{*arguments};

        auto hash_view = [&] ()
        {
            if constexpr (!std::ranges::view<view_t>)
            {
                return seqan3::views::minimiser_hash(arguments->shape,
                                                     seqan3::window_size{arguments->window_size},
                                                     seqan3::seed{adjust_seed(arguments->shape_weight)});
            }
            else
            {
                return seqan3::views::minimiser_hash(arguments->shape,
                                                     seqan3::window_size{arguments->window_size},
                                                     seqan3::seed{adjust_seed(arguments->shape_weight)})
                       | hash_filter_view;
            }
        };

        if (arguments->bin_path.size() == 1)
        {
            metadata meta(arguments->ref_meta_path);

            auto & ibf = index.ibf();
            size_t i = 0;
            for (auto && [seq] : sequence_file_t{arguments->bin_path[0][0]})
            {
                for (auto & seg : meta.segments_from_ind(i))
                {
                    for (auto && value : seq | seqan3::views::slice(seg.start, seg.start + seg.len) | hash_view())
                        ibf.emplace(value, seqan3::bin_index{seg.id});
                }
                i++;
            }
        }
        else
        {
            auto clustered_reference_worker = [&] (auto && zipped_view, auto &&)
            {
                auto & ibf = index.ibf();

                for (auto && [file_names, bin_number] : zipped_view)
                    for (auto && file_name : file_names)
                        for (auto && [seq] : sequence_file_t{file_name})
                            for (auto && value : seq | hash_view())
                                ibf.emplace(value, seqan3::bin_index{bin_number});
            };
            call_parallel_on_bins(clustered_reference_worker, arguments->bin_path, arguments->threads);
        }
        return index;
    }
};

} // namespace valik
