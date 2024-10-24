// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements seqan::hibf::bin_size_in_bits.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <utilities/prepare/compute_bin_size.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace raptor
{

std::filesystem::path get_output_path(std::filesystem::path const & output_dir, std::filesystem::path const & file_name)
{
    std::filesystem::path result{output_dir};
    result /= file_name.stem();
    result += ".dummy_extension"; // https://github.com/seqan/raptor/issues/355
    return result;
}

std::filesystem::path get_output_path(std::filesystem::path const & output_dir, std::filesystem::path const & file_name, size_t const bin_id)
{
    std::filesystem::path result{output_dir};
    result /= file_name.stem();
    result += ".";
    result += std::to_string(bin_id);
    result += ".dummy_extension"; // https://github.com/seqan/raptor/issues/355
    return result;
}

void compute_minimiser(valik::build_arguments const & arguments)
{
    if (arguments.bin_path.size() > 1)
    {
        file_reader<file_types::sequence> const reader{arguments.shape, arguments.window_size};
        auto cluster_worker = [&](auto && zipped_view, auto &&)
        {
            for (auto && [file_names, bin_number] : zipped_view)
            {
                std::filesystem::path const file_name{file_names[0]};
                size_t const seq_size = std::filesystem::file_size(file_name);
                std::filesystem::path output_path = get_output_path(arguments.out_dir, file_name);

                std::filesystem::path const minimiser_file =
                    std::filesystem::path{output_path}.replace_extension("minimiser");
                std::filesystem::path const progress_file =
                    std::filesystem::path{output_path}.replace_extension("in_progress");
                std::filesystem::path const header_file = std::filesystem::path{output_path}.replace_extension("header");
            
                // If we are already done with this file, we can skip it. Otherwise, we create a ".in_progress" file to keep
                // track of whether the minimiser computation was successful.
                bool const already_done = std::filesystem::exists(header_file) && std::filesystem::exists(header_file) && 
                                          !std::filesystem::exists(progress_file);

                if (already_done)
                    continue;
                else
                    std::ofstream outfile{progress_file, std::ios::binary};

                std::unordered_set<uint64_t> distinct_minimisers{};
                // The set is (re-)constructed for each file. The alternative is to construct it once for each thread
                // and clear+reuse it for every file that a thread works on. However, this dramatically increases
                // memory consumption because the map will stay as big as needed for the biggest encountered file.

                reader.for_each_hash(file_names,
                                     [&](auto && hash)
                                     {
                                         distinct_minimisers.insert(hash);
                                     });

                uint64_t count{};

                {
                    //!TODO: apply k-mer count cutoffs in metagenome search
                    std::ofstream outfile{minimiser_file, std::ios::binary};
                    for (auto && hash : distinct_minimisers)
                    {
                        outfile.write(reinterpret_cast<const char *>(&hash), sizeof(hash));
                        ++count;
                    }
                }

                {
                    std::ofstream headerfile{header_file};
                    headerfile << arguments.shape.to_string() << '\t' << arguments.window_size << '\t' << count << '\t' << seq_size << '\n';
                }

                std::filesystem::remove(progress_file);
            }
        };
        
        size_t const chunk_size =
            std::max<size_t>(1, std::floor(arguments.bin_path.size() / static_cast<double>(arguments.threads)));
        auto chunked_view = seqan3::views::zip(arguments.bin_path, std::views::iota(0u)) | seqan3::views::chunk(chunk_size);
        seqan3::detail::execution_handler_parallel executioner{arguments.threads};
        executioner.bulk_execute(std::move(cluster_worker), std::move(chunked_view), []() {});
    }
    else
    {
        valik::metadata meta(arguments.ref_meta_path);
        std::filesystem::path const file_name{arguments.bin_path[0][0]};
        
        auto segment_worker = [&](const auto && zipped_view, auto &&)
        {
            for (auto && [shared_record, bin_number] : zipped_view)
            {
                //!TODO: bin number is not used. How to construct view of shared records to fulfill execution handler requirements? 
                //seqan3::debug_stream << bin_number << '\n';
                auto const & seg = shared_record.segment;

                std::filesystem::path output_path = get_output_path(arguments.out_dir, file_name, seg.id);
                std::filesystem::path const minimiser_file =
                std::filesystem::path{output_path}.replace_extension("minimiser");
                std::filesystem::path const progress_file =
                    std::filesystem::path{output_path}.replace_extension("in_progress");
                std::filesystem::path const header_file = std::filesystem::path{output_path}.replace_extension("header");

                // If we are already done with this file, we can skip it. Otherwise, we create a ".in_progress" file to keep
                // track of whether the minimiser computation was successful.
                bool const already_done = std::filesystem::exists(header_file) && std::filesystem::exists(header_file) && 
                                          !std::filesystem::exists(progress_file);

                if (already_done)
                    continue;
                else
                    std::ofstream outfile{progress_file, std::ios::binary};

                std::unordered_map<uint64_t, uint8_t> minimiser_table{};
                // The map is (re-)constructed for each segment. The alternative is to construct it once for each thread
                // and clear+reuse it for every file that a thread works on. However, this dramatically increases
                // memory consumption because the map will stay as big as needed for the biggest encountered file.

                auto hash_view = [&] ()
                {
                    return seqan3::views::minimiser_hash(arguments.shape,
                                                         seqan3::window_size{arguments.window_size},
                                                         seqan3::seed{adjust_seed(arguments.shape_weight)});
                };

                for (auto const value : *shared_record.underlying_sequence | seqan3::views::slice(seg.start, seg.start + seg.len) | hash_view())
                {
                    minimiser_table[value] = std::min<uint8_t>(254u, minimiser_table[value] + 1);  
                }

                uint64_t count{};
                {
                    std::ofstream outfile{minimiser_file, std::ios::binary};
                    for (auto && [hash, occurrences] : minimiser_table)
                    {
                        if (occurrences >= arguments.kmer_count_min_cutoff && occurrences <= arguments.kmer_count_max_cutoff)
                        {
                            outfile.write(reinterpret_cast<const char *>(&hash), sizeof(hash));
                            ++count;
                        }
                    }
                }

                {
                    std::ofstream headerfile{header_file};
                    headerfile << arguments.shape.to_string() << '\t' << arguments.window_size << '\t' << count << '\t' << seg.len << '\n';
                }

                std::filesystem::remove(progress_file);
            }
        };

        using sequence_t = seqan3::dna4_vector;
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
        std::vector<valik::shared_reference_record<sequence_t>> reference_records;
        reference_records.reserve(arguments.threads);
        size_t seq_ind{0};
        size_t processed_bin_count{0};
        for (auto && [seq] : sequence_file_t{arguments.bin_path[0][0]})
        {
            auto shared_seq = std::make_shared<sequence_t>(seq);
            for (auto && seg : meta.segments_from_ind(seq_ind))
            {
                reference_records.emplace_back(shared_seq, std::move(seg));
                if (reference_records.size() == arguments.threads)
                {
                    auto chunked_view = seqan3::views::zip(reference_records, std::views::iota(processed_bin_count)) | seqan3::views::chunk(1u);
                    seqan3::detail::execution_handler_parallel executioner{arguments.threads};
                    executioner.bulk_execute(std::move(segment_worker), std::move(chunked_view), []() {});
                    processed_bin_count += reference_records.size();
                    reference_records.clear();
                }  
            }
            seq_ind++;
        }

        if (!reference_records.empty())
        {
            seqan3::detail::execution_handler_parallel executioner{reference_records.size()};
            auto chunked_view = seqan3::views::zip(reference_records, std::views::iota(processed_bin_count)) | seqan3::views::chunk(1u);
            executioner.bulk_execute(std::move(segment_worker), std::move(chunked_view), []() {});
        }
    }
}

namespace detail
{

size_t kmer_count_from_minimiser_files(std::vector<std::vector<std::string>> const & minimiser_bin_path, uint8_t const threads)
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

    valik::call_parallel_on_bins(worker, minimiser_bin_path, threads);

    std::string shape_string{};
    uint64_t window_size{};
    size_t max_count{};
    uint64_t bin_size{};

    biggest_file.replace_extension("header");
    std::ifstream file_stream{biggest_file};
    file_stream >> shape_string >> window_size >> max_count >> bin_size;

    return max_count;
}

size_t kmer_count_from_sequence_files(valik::build_arguments const & arguments)
{
    size_t max_count{};
    std::mutex callback_mutex{};

    auto callback = [&callback_mutex, &max_count](size_t const count)
    {
        std::lock_guard<std::mutex> guard{callback_mutex};
        max_count = std::max<size_t>(max_count, count);
    };
    
    if (arguments.bin_path.size() > 1)
    {
        file_reader<file_types::sequence> const reader{arguments.shape, arguments.window_size};

        auto cluster_worker = [&callback, &reader](auto && zipped_view, auto &&)
        {
            std::unordered_set<uint64_t> kmers;
            auto insert_it = std::inserter(kmers, kmers.end());
            for (auto && [file_names, bin_number] : zipped_view)
            {
                kmers.clear();
                reader.hash_into(file_names, insert_it);
                callback(kmers.size());
            }   
        };

        // callback max_count which is the max number of k-mers in any sequence 
        valik::call_parallel_on_bins(cluster_worker, arguments.bin_path, arguments.threads);
    }
    else
    {
        using sequence_t = seqan3::dna4_vector;
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
        std::vector<valik::shared_reference_record<sequence_t>> reference_records;
        reference_records.reserve(arguments.threads);
        size_t seq_ind{0};
        size_t processed_bin_count{0};

        auto hash_view = [&] ()
        {
            return seqan3::views::minimiser_hash(arguments.shape,
                                                 seqan3::window_size{arguments.window_size},
                                                 seqan3::seed{adjust_seed(arguments.shape.count())});
        };

        auto segment_worker = [&callback, &hash_view](const auto && zipped_view, auto &&)
        {
            std::unordered_set<uint64_t> kmers;
            auto insert_it = std::inserter(kmers, kmers.end());
            for (auto && [shared_record, bin_number] : zipped_view)
            {
                kmers.clear();
                //!TODO: bin number is not used. How to construct view of shared records to fulfill execution handler requirements? 
                //seqan3::debug_stream << bin_number << '\n';
                auto const & seg = shared_record.segment;
                std::ranges::copy(*shared_record.underlying_sequence | 
                                  seqan3::views::slice(seg.start, seg.start + seg.len) | hash_view(), insert_it);

                callback(kmers.size());
            }
        };

        for (auto && [seq] : sequence_file_t{arguments.bin_path[0][0]})
        {
            valik::metadata meta(arguments.ref_meta_path);
            auto shared_seq = std::make_shared<sequence_t>(seq);
            for (auto && seg : meta.segments_from_ind(seq_ind))
            {
                reference_records.emplace_back(shared_seq, std::move(seg));
                if (reference_records.size() == arguments.threads)
                {
                    auto chunked_view = seqan3::views::zip(reference_records, std::views::iota(processed_bin_count)) | seqan3::views::chunk(1u);
                    seqan3::detail::execution_handler_parallel executioner{arguments.threads};
                    executioner.bulk_execute(std::move(segment_worker), std::move(chunked_view), []() {});
                    processed_bin_count += reference_records.size();
                    reference_records.clear();
                }  
            }
            seq_ind++;
        }

        if (!reference_records.empty())
        {
            seqan3::detail::execution_handler_parallel executioner{reference_records.size()};
            auto chunked_view = seqan3::views::zip(reference_records, std::views::iota(processed_bin_count)) | seqan3::views::chunk(1u);
            executioner.bulk_execute(std::move(segment_worker), std::move(chunked_view), []() {});
        }
    }

    return max_count;
}

} // namespace detail


size_t compute_bin_size(valik::build_arguments const & arguments)
{
    std::vector<std::vector<std::string>> minimiser_files{};
    if (arguments.input_is_minimiser)
    {
        minimiser_files = parse_bin_paths(arguments);
    }

    size_t max_count = arguments.input_is_minimiser
                               ? detail::kmer_count_from_minimiser_files(minimiser_files, arguments.threads)
                               : detail::kmer_count_from_sequence_files(arguments);

    assert(max_count > 0u);

    return seqan::hibf::build::bin_size_in_bits(
        {.fpr = arguments.fpr, .hash_count = arguments.hash, .elements = max_count});
}

} // namespace raptor
