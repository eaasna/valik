#pragma once

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include <valik/shared.hpp>

namespace valik
{

template <typename algorithm_t>
void call_parallel_on_bins(algorithm_t && worker,
                           std::vector<std::vector<std::string>> const & bin_paths,
                           uint8_t const threads)
{
    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(bin_paths.size() / threads),
                                                 8u,
                                                 64u); // each thread will process between 8 to 64 bins
    auto chunked_view = seqan3::views::zip(bin_paths, std::views::iota(0u)) | // pair each database input file in the bin_path vector with bin number
                        seqan3::views::chunk(chunk_size); // how many input files are processed by each thread
    seqan3::detail::execution_handler_parallel executioner{threads};
    executioner.bulk_execute(std::move(worker), std::move(chunked_view), [](){});
}

} // namespace valik
