// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements seqan::hibf::bin_size_in_bits.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <utilities/prepare/compute_bin_size.hpp>

namespace raptor
{

std::filesystem::path get_output_path(std::filesystem::path const & output_dir, std::filesystem::path const & file_name)
{
    std::filesystem::path result{output_dir};
    bool const is_compressed = raptor::cutoff::file_is_compressed(file_name);
    result /= is_compressed ? file_name.stem().stem() : file_name.stem();
    result += ".dummy_extension"; // https://github.com/seqan/raptor/issues/355
    return result;
}

void compute_minimiser(valik::build_arguments const & arguments)
{
    file_reader<file_types::sequence> const reader{arguments.shape, arguments.window_size};
    raptor::cutoff const cutoffs{arguments};

    auto worker = [&](auto && zipped_view, auto &&)
    {
        for (auto && [file_names, bin_number] : zipped_view)
        {
            std::filesystem::path const file_name{file_names[0]};
            std::filesystem::path output_path = get_output_path(arguments.out_dir, file_name);

            std::filesystem::path const minimiser_file =
                std::filesystem::path{output_path}.replace_extension("minimiser");
            std::filesystem::path const progress_file =
                std::filesystem::path{output_path}.replace_extension("in_progress");
            std::filesystem::path const header_file = std::filesystem::path{output_path}.replace_extension("header");
            
            // If we are already done with this file, we can skip it. Otherwise, we create a ".in_progress" file to keep
            // track of whether the minimiser computation was successful.
            bool const already_done = std::filesystem::exists(header_file) && !std::filesystem::exists(progress_file);

            if (already_done)
                continue;
            else
                std::ofstream outfile{progress_file, std::ios::binary};

            // The hash table stores how often a minimiser appears. It does not matter whether a minimiser appears
            // 50 times or 2000 times, it is stored regardless because the biggest cutoff value is 50. Hence,
            // the hash table stores only values up to 254 to save memory.
            std::unordered_map<uint64_t, uint8_t> minimiser_table{};
            // The map is (re-)constructed for each file. The alternative is to construct it once for each thread
            // and clear+reuse it for every file that a thread works on. However, this dramatically increases
            // memory consumption because the map will stay as big as needed for the biggest encountered file.

            reader.for_each_hash(file_names,
                                 [&](auto && hash)
                                 {
                                     minimiser_table[hash] = std::min<uint8_t>(254u, minimiser_table[hash] + 1);
                                 });

            uint8_t const cutoff = cutoffs.get(file_name);
            uint64_t count{};

            {
                std::ofstream outfile{minimiser_file, std::ios::binary};
                for (auto && [hash, occurrences] : minimiser_table)
                {
                    if (occurrences >= cutoff)
                    {
                        outfile.write(reinterpret_cast<const char *>(&hash), sizeof(hash));
                        ++count;
                    }
                }
            }

            {
                std::ofstream headerfile{header_file};
                headerfile << arguments.shape.to_string() << '\t' << arguments.window_size << '\t'
                           << static_cast<uint16_t>(cutoff) << '\t' << count << '\n';
            }

            std::filesystem::remove(progress_file);
        }
    };

    size_t const chunk_size =
        std::max<size_t>(1, std::floor(arguments.bin_path.size() / static_cast<double>(arguments.threads)));
    auto chunked_view = seqan3::views::zip(arguments.bin_path, std::views::iota(0u)) | seqan3::views::chunk(chunk_size);
    seqan3::detail::execution_handler_parallel executioner{arguments.threads};
    executioner.bulk_execute(std::move(worker), std::move(chunked_view), []() {});
}

namespace detail
{

size_t kmer_count_from_minimiser_files(std::vector<std::vector<std::string>> const & bin_path, uint8_t const threads)
{
    std::mutex callback_mutex{};
    size_t max_filesize{};
    std::filesystem::path biggest_file{};

    auto callback = [&callback_mutex, &biggest_file, &max_filesize](std::filesystem::path const path, size_t const size)
    {
        std::lock_guard<std::mutex> guard{callback_mutex};
        if (size > max_filesize)
        {
            max_filesize = size;
            biggest_file = path;
        }
    };

    auto worker = [&callback](auto && zipped_view, auto &&)
    {
        std::filesystem::path minimiser_file{};
        std::filesystem::path biggest_file{};
        size_t max_filesize{};

        for (auto && [file_names, bin_number] : zipped_view)
        {
            for (auto && file_name : file_names)
            {
                minimiser_file = file_name;
                size_t const size = std::filesystem::file_size(minimiser_file);
                if (size > max_filesize)
                {
                    max_filesize = size;
                    biggest_file = minimiser_file;
                }
            }
        }

        callback(biggest_file, max_filesize);
    };

    valik::call_parallel_on_bins(worker, bin_path, threads);

    std::string shape_string{};
    uint64_t window_size{};
    uint16_t cutoff{};
    size_t max_count{};

    biggest_file.replace_extension("header");
    std::ifstream file_stream{biggest_file};
    file_stream >> shape_string >> window_size >> cutoff >> max_count;

    return max_count;
}

size_t kmer_count_from_sequence_files(std::vector<std::vector<std::string>> const & bin_path,
                                      uint8_t const threads,
                                      seqan3::shape const & shape,
                                      uint32_t const window_size)
{
    size_t max_count{};
    std::mutex callback_mutex{};
    file_reader<file_types::sequence> const reader{shape, window_size};

    auto callback = [&callback_mutex, &max_count](auto && view)
    {
        auto const count = std::ranges::distance(view);
        {
            std::lock_guard<std::mutex> guard{callback_mutex};
            max_count = std::max<size_t>(max_count, count);
        }
    };

    auto worker = [&callback, &reader](auto && zipped_view, auto &&)
    {
        for (auto && [file_names, bin_number] : zipped_view)
            reader.on_hash(file_names, callback);
    };

    // callback max_count which is the max number of k-mers in any sequence 
    valik::call_parallel_on_bins(worker, bin_path, threads);

    return max_count;
}

} // namespace detail


size_t compute_bin_size(valik::build_arguments const & arguments)
{
    size_t max_count = arguments.input_is_minimiser
                               ? detail::kmer_count_from_minimiser_files(arguments.bin_path, arguments.threads)
                               : detail::kmer_count_from_sequence_files(arguments.bin_path,
                                                                        arguments.threads,
                                                                        arguments.shape,
                                                                        arguments.window_size);

    if (arguments.bin_path.size() == 1)
        max_count /= arguments.bins;

    assert(max_count > 0u);

    return seqan::hibf::build::bin_size_in_bits(
        {.fpr = arguments.fpr, .hash_count = arguments.hash, .elements = max_count});
}

} // namespace raptor
