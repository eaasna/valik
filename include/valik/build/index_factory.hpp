#pragma once

#include <utilities/prepare/parse_bin_paths.hpp>
#include <valik/build/call_parallel_on_bins.hpp>
#include <valik/index.hpp>
#include <valik/split/metadata.hpp>

#include <seqan3/search/views/minimiser_hash.hpp>

namespace valik
{

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
    [[nodiscard]] auto operator()() const
    {
        return construct();
    }

private:
    build_arguments const * const arguments{nullptr};

    auto construct() const
    {
        assert(arguments != nullptr);

        valik_index<> index{*arguments};
        auto & ibf = index.ibf();
        auto & entropy_ranking = index.entropy_ranking();
        std::vector<std::pair<size_t, double>> entropy_map{};

        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
        auto hash_view = [&] ()
        {
            return seqan3::views::minimiser_hash(arguments->shape,
                                                 seqan3::window_size{arguments->window_size},
                                                 seqan3::seed{adjust_seed(arguments->shape_weight)});
        };

        if (arguments->input_is_minimiser)
        {
            auto minimiser_worker = [&] (auto && zipped_view, auto &&)
            {
                for (auto && [file_name, bin_number] : zipped_view)
                {
                    std::ifstream fin{file_name, std::ios::binary};
                    uint64_t value;
                    while (fin.read(reinterpret_cast<char *>(&value), sizeof(value)))
                    {
                        ibf.emplace(value, seqan3::bin_index{bin_number});
                    }
                }
            };

            std::vector<std::string> file_paths = parse_bin_paths(*arguments);
            call_parallel_on_bins(minimiser_worker, file_paths, arguments->threads);

            std::vector<std::string> header_paths = parse_bin_paths(*arguments, "header");
            std::string shape_string{};
            uint64_t window_size{};
            size_t count{};
            uint64_t bin_size{};
            entropy_ranking.reserve(header_paths.size());
            for (auto && [file_name, bin_number] : seqan3::views::zip(header_paths, std::views::iota(0u)))
            {
                std::ifstream file_stream{file_name};
                file_stream >> shape_string >> window_size >> count >> bin_size;
                entropy_map.emplace_back(std::make_pair((size_t) bin_number, (double) count / (double) bin_size));
            }

            std::ranges::sort(entropy_map.begin(), entropy_map.end(), [](const std::pair<size_t, double> &a, const std::pair<size_t, double> &b)
            { 
                return a.second > b.second; 
            });

            for (auto && [bin_id, entropy] : entropy_map)
                entropy_ranking.push_back(bin_id);
        }
        else if (arguments->bin_path.size() > 1)
        {
            auto clustered_reference_worker = [&] (auto && zipped_view, auto &&)
            {
                for (auto && [file_name, bin_number] : zipped_view)
                {
                    for (auto && record : sequence_file_t{file_name})
                    {
                        for (auto && value : record.sequence() | hash_view())
                        {
                            ibf.emplace(value, seqan3::bin_index{bin_number});
                        }
                    }
                }
            };
            
            call_parallel_on_bins(clustered_reference_worker, arguments->bin_path, arguments->threads);
        }
        else
        {
            // if no .minimiser files exist then a single sequence file will be scanned sequentially to build the IBF
            metadata meta(arguments->ref_meta_path);

            auto min_view = hash_view();

            size_t i{0};
            for (auto && [seq] : sequence_file_t{arguments->bin_path[0]})
            {
                for (auto & seg : meta.segments_from_ind(i))
                {
                    for (auto && value : seq | seqan3::views::slice(seg.start, seg.start + seg.len) | min_view)
                        ibf.emplace(value, seqan3::bin_index{seg.id});
                }
                i++;
            }
        }

        return index;
    }
};

} // namespace valik
