#pragma once

#include <valik/search/cart_query_io.hpp>
#include <valik/search/iterate_queries.hpp>
#include <valik/search/load_index.hpp>
#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>
#include <utilities/consolidate/merge_processes.hpp>
#include <utilities/cart_queue.hpp>

namespace valik::app
{

/**
 * @brief Function that calls Valik prefiltering and launches parallel processes of Stellar search.
 *
 * @tparam compressed IBF layout type.
 * @param arguments Command line arguments.
 * @param time_statistics Run-time statistics.
 * @return false if search failed.
 */
template <bool compressed>
bool search_distributed(search_arguments const & arguments, search_time_statistics & time_statistics)
{
    using index_structure_t = std::conditional_t<compressed, index_structure::ibf_compressed, index_structure::ibf>;
    auto index = valik_index<index_structure_t>{};

    {
        auto start = std::chrono::high_resolution_clock::now();
        load_index(index, arguments.index_file);
        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    }

    std::optional<metadata> ref_meta;
    if (!arguments.ref_meta_path.empty())
        ref_meta = metadata(arguments.ref_meta_path);

    env_var_pack var_pack{};
    auto queue = cart_queue<query_record>{index.ibf().bin_count(), arguments.cart_max_capacity, arguments.max_queued_carts};

    std::mutex mutex;
    execution_metadata exec_meta(arguments.threads);

    bool error_in_search = false; // indicates if an error happen inside this lambda
    auto consumerThreads = std::vector<std::jthread>{};
    for (size_t threadNbr = 0; threadNbr < arguments.threads; ++threadNbr)
    {
        consumerThreads.emplace_back(
        [&, threadNbr]() {
            auto& thread_meta = exec_meta.table[threadNbr];
            // this will block until producer threads have added carts to queue
            for (auto next = queue.dequeue(); next; next = queue.dequeue())
            {
                auto & [bin_id, records] = *next;

                std::unique_lock g(mutex);
                std::filesystem::path cart_queries_path = var_pack.tmp_path / std::string("query_" + std::to_string(bin_id) +
                                                          "_" + std::to_string(exec_meta.bin_count[bin_id]++) + ".fasta");
                g.unlock();

                thread_meta.output_files.push_back(cart_queries_path.string() + ".gff");

                write_cart_queries(records, cart_queries_path);

                std::vector<std::string> process_args{};
                process_args.insert(process_args.end(), {var_pack.stellar_exec, "--version-check", "0", "--verbose", "-a", "dna"});

                if (ref_meta)
                {
                    // search segments of a single reference file
                    auto ref_len = ref_meta->total_len;
                    auto seg = ref_meta->segment_from_bin(bin_id);
                    process_args.insert(process_args.end(), {index.bin_path()[0][0], std::string(cart_queries_path),
                                                            "--referenceLength", std::to_string(ref_len),
                                                            "--sequenceOfInterest", std::to_string(seg.seq_vec[0]),
                                                            "--segmentBegin", std::to_string(seg.start),
                                                            "--segmentEnd", std::to_string(seg.start + seg.len)});
                }
                else
                {
                    // search a reference database of bin sequence files
                    if (index.bin_path().size() < (size_t) bin_id) {
                        throw std::runtime_error("Could not find reference file with index " + std::to_string(bin_id) +
                        ". Did you forget to provide metadata to search segments in a single reference file instead?");
                    }
                    process_args.insert(process_args.end(), {index.bin_path()[bin_id][0], std::string(cart_queries_path)});
                }

                if (arguments.write_time)
                    process_args.insert(process_args.end(), "--time");

                // ==========================================
                //!WORKAROUND: Stellar does not allow smaller error rates
                // ==========================================
                float numEpsilon = std::max(arguments.error_rate, (float) 0.00001);
                process_args.insert(process_args.end(), {"-e", std::to_string(numEpsilon),
                                                        "-l", std::to_string(arguments.pattern_size),
                                                        "-o", std::string(cart_queries_path) + ".gff"});
                
                process_args.insert(process_args.end(), {"--repeatPeriod", std::to_string(arguments.maxRepeatPeriod)});
                process_args.insert(process_args.end(), {"--repeatLength", std::to_string(arguments.minRepeatLength)});
                process_args.insert(process_args.end(), {"--verification", arguments.strVerificationMethod});
                process_args.insert(process_args.end(), {"--xDrop", std::to_string(arguments.xDrop)});
                process_args.insert(process_args.end(), {"--abundanceCut", std::to_string(arguments.qgramAbundanceCut)});

                auto start = std::chrono::high_resolution_clock::now();
                external_process process(process_args);
                auto end = std::chrono::high_resolution_clock::now();

                thread_meta.time_statistics.emplace_back(0.0 + std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count());

                thread_meta.text_out << process.cout();
                thread_meta.text_out << process.cerr();

                if (process.status() != 0) {
                    std::unique_lock g(mutex); // make sure that our output is synchronized
                    std::cerr << "error running VALIK_STELLAR\n";
                    std::cerr << "call:";
                    for (auto args : process_args) {
                        std::cerr << " " << args;
                    }
                    std::cerr << '\n';
                    std::cerr << process.cerr() << '\n';
                    error_in_search = true;
                }
            }
        });
    }

    auto start = std::chrono::high_resolution_clock::now();
    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};
    iterate_distributed_queries(arguments, index.ibf(), thresholder, queue);

    queue.finish(); // Flush carts that are not empty yet
    consumerThreads.clear();
    auto end = std::chrono::high_resolution_clock::now();
    time_statistics.search_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    // merge output files and metadata from threads
    bool error_in_merge = merge_processes(arguments, time_statistics, exec_meta, var_pack);

    return error_in_search && error_in_merge;
}

} // namespace valik::app
