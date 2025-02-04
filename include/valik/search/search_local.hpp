#pragma once

#include "future"

#include <valik/search/cart_query_io.hpp>
#include <valik/search/iterate_queries.hpp>
#include <valik/search/load_index.hpp>
#include <valik/shared.hpp>
#include <utilities/cart_queue.hpp>
#include <utilities/consolidate/merge_processes.hpp>
#include <utilities/threshold/search_kmer_profile.hpp>
#include <utilities/threshold/filtering_request.hpp>

#include <dream_stellar/database_id_map.hpp>
#include <dream_stellar/diagnostics/print.tpp>
#include <dream_stellar/stellar_launcher.hpp>
#include <dream_stellar/stellar_output.hpp>
#include <dream_stellar/io/import_sequence.hpp>
#include <dream_stellar/utils/stellar_app_runtime.hpp>
#include <dream_stellar/shared.hpp>

namespace valik::app
{

template <typename TSize, typename id_t>
static inline dream_stellar::StellarOptions make_thread_options(search_arguments const & arguments, metadata const & ref_meta, 
                                                          std::filesystem::path const & cart_queries_path, TSize const refLen, id_t const bin_id)
{
    dream_stellar::StellarOptions threadOptions{};
    threadOptions.alphabet = "dna";            // one of: dna, rna, protein, char
    threadOptions.verbose = arguments.very_verbose;
    threadOptions.queryFile = "in memory";
    threadOptions.prefilteredSearch = true;
    threadOptions.referenceLength = refLen;
    threadOptions.searchSegment = true;
    auto seg = ref_meta.segment_from_bin(bin_id);
    threadOptions.binSequences.emplace_back(seg.seq_vec[0]);
    threadOptions.segmentBegin = seg.start;
    threadOptions.segmentEnd = seg.start + seg.len;
    threadOptions.minLength = arguments.minLength;
    threadOptions.epsilon = dream_stellar::utils::fraction::from_double_with_limit(arguments.error_rate, arguments.minLength).limit_denominator();
    threadOptions.outputFile = cart_queries_path.string() + ".gff";
    
    {
        threadOptions.maxRepeatPeriod = arguments.maxRepeatPeriod;
        threadOptions.minRepeatLength = arguments.minRepeatLength;
        threadOptions.strVerificationMethod = arguments.strVerificationMethod;
        threadOptions.xDrop = arguments.xDrop;
        threadOptions.qgramAbundanceCut = arguments.qgramAbundanceCut;
        threadOptions.numMatches = arguments.numMatches;
        threadOptions.compactThresh = arguments.compactThresh;
    }

    return threadOptions;
}

/**
 * @brief Function that calls Valik prefiltering and launches parallel threads of Stellar search.
 *
 * @tparam is_split Split query sequences.
 * @param arguments Command line arguments.
 * @param time_statistics Run-time statistics.
 * @return false if search failed.
 */
template <bool is_split, bool stellar_only>
bool search_local(search_arguments & arguments, search_time_statistics & time_statistics)
{
    using index_structure_t = index_structure::ibf;
    auto index = valik_index<index_structure_t>{};

    if (!stellar_only)
    {
        auto start = std::chrono::high_resolution_clock::now();
        load_index(index, arguments.index_file);
        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();        
    }

    metadata ref_meta = metadata(arguments.ref_meta_path);
    if (stellar_only)
    {
        arguments.bin_path.clear(); // in case inserted from index
        for (auto & f : ref_meta.files)
            arguments.bin_path.push_back(std::vector<std::string>{f.path});

        auto prefilter_bin_count = ref_meta.seg_count;
        split_arguments stellar_dist_arguments;
        // stellar search without prefiltering
        
        // for some number of reference sequences split sequences into as many segments as is the next multiple of thread count
        if (ref_meta.seq_count % arguments.threads > 0)
            stellar_dist_arguments.seg_count = ref_meta.seq_count + (arguments.threads - ref_meta.seq_count % arguments.threads);
        else
            stellar_dist_arguments.seg_count = ref_meta.seq_count;

        stellar_dist_arguments.pattern_size = ref_meta.pattern_size;
        ref_meta.update_segments_for_distributed_stellar(stellar_dist_arguments);
        // update cart queue parameters for distributed stellar
        arguments.cart_max_capacity = std::round(arguments.cart_max_capacity * prefilter_bin_count / (double) ref_meta.seg_count);
        arguments.max_queued_carts = std::round(arguments.max_queued_carts / (double) prefilter_bin_count * ref_meta.seg_count);
    }

    if (arguments.max_queued_carts == std::numeric_limits<uint32_t>::max()) // if no user input
        arguments.max_queued_carts = ref_meta.seg_count;

    env_var_pack var_pack{};
    std::optional<metadata> query_meta;
    if (arguments.split_query && !stellar_only)
    {
        
        query_meta = metadata(arguments);
        if (arguments.verbose)
        {
            std::cout << "\n-----------Preprocessing queries-----------\n";
            std::cout << "database size " << query_meta.value().total_len << "bp\n";
            std::cout << "segment count " << query_meta.value().seg_count << '\n';
            std::cout << "segment len " << std::to_string((uint64_t) std::round(query_meta.value().total_len / 
                                                         (double) query_meta.value().seg_count)) << "bp\n";
        }

        if (!arguments.manual_parameters)
        {
            search_pattern pattern(arguments.errors, arguments.pattern_size);
            param_space space;
            param_set params(arguments.shape, arguments.threshold);
            filtering_request request(pattern, ref_meta, query_meta.value());
            if ((request.fpr(params) > 0.2) && (arguments.search_type != search_kind::STELLAR))
                std::cerr << "WARNING: Prefiltering will be inefficient for a high error rate.\n";

            if (arguments.verbose)
            {
                std::cout.precision(3);

                std::cout << "\n-----------Search parameters-----------\n";
                if (arguments.shape_size == arguments.shape_weight)
                    std::cout << "kmer size " << std::to_string(arguments.shape_size) << '\n';
                else 
                    std::cout << "kmer shape " << arguments.shape.to_string() << '\n';
                std::cout << "window size " << std::to_string(arguments.window_size) << '\n';
                switch (arguments.search_type)
                {
                    case search_kind::LEMMA: std::cout << "k-mer lemma "; break;
                    case search_kind::MINIMISER: std::cout << "minimiser "; break;
                    case search_kind::HEURISTIC: std::cout << "heuristic "; break;
                    default: break;
                }
                std::cout << "threshold ";
                std::cout << std::to_string(arguments.threshold) << '\n';

                std::cout << "FNR " << arguments.fnr << '\n';
                std::cout << "FPR " << request.fpr(params) << '\n';
            }
        }
    }

    using TAlphabet = seqan2::Dna;
    using TSequence = seqan2::String<TAlphabet>;
    // the queue hands records over from the producer threads (valik prefiltering) to the consumer threads (stellar search) 
    auto queue = cart_queue<shared_query_record<TSequence>>{ref_meta.seg_count, arguments.cart_max_capacity, arguments.max_queued_carts};

    std::mutex mutex;
    execution_metadata exec_meta(arguments.threads);

    // negative (reverse complemented) database strand
    bool const reverse = true /*threadOptions.reverse && threadOptions.alphabet != "protein" && threadOptions.alphabet != "char" */;
    seqan2::StringSet<TSequence> databases;
    using TSize = decltype(length(databases[0]));
    seqan2::StringSet<TSequence> reverseDatabases;
    seqan2::StringSet<seqan2::CharString> databaseIDs;
    std::ofstream disabledQueriesFile;
    TSize refLen;

    if (arguments.disableThresh != std::numeric_limits<size_t>::max())
    {
        std::filesystem::path disabledPath = arguments.out_file;
        disabledPath.replace_extension(".disabled.fa");
        disabledQueriesFile.open(((std::string) disabledPath).c_str(), ::std::ios_base::out | ::std::ios_base::app);
        if (!disabledQueriesFile.is_open())
        {
            std::cerr << "Could not open file for disabled queries." << std::endl;
            return false;
        }
    }

    dream_stellar::stellar_runtime input_databases_time{};
    bool const databasesSuccess = input_databases_time.measure_time([&]()
    {
        std::cout << "Launching stellar search on a shared memory machine...\n";
        return dream_stellar::_importAllSequences(arguments.bin_path[0][0].c_str(), "database", databases, databaseIDs, refLen, std::cout, std::cerr);
    });
    if (!databasesSuccess)
        return false;

    if (reverse)
    {
        for (auto database : databases)
        {
            reverseComplement(database);
            seqan2::appendValue(reverseDatabases, database, seqan2::Generous());
        }
    }

    time_statistics.ref_io_time += input_databases_time.milliseconds() / 1000;
    dream_stellar::DatabaseIDMap<TAlphabet> databaseIDMap{databases, databaseIDs};
    dream_stellar::DatabaseIDMap<TAlphabet> reverseDatabaseIDMap{reverseDatabases, databaseIDs};

    bool error_in_search = false; // indicates if an error happened inside this lambda
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

                dream_stellar::stellar_app_runtime stellarThreadTime{};
                auto current_time = stellarThreadTime.now();
                dream_stellar::StellarOptions threadOptions = make_thread_options(arguments, ref_meta, cart_queries_path, refLen, bin_id);

                using TDatabaseSegment = dream_stellar::StellarDatabaseSegment<TAlphabet>;
                using TQuerySegment = seqan2::Segment<seqan2::String<TAlphabet> const, seqan2::InfixSegment>;

                dream_stellar::_writeFileNames(threadOptions, thread_meta.text_out);
                dream_stellar::_writeSpecifiedParams(threadOptions, thread_meta.text_out);
                dream_stellar::_writeCalculatedParams(threadOptions, thread_meta.text_out);   // calculate qGram
                thread_meta.text_out << std::endl;

                // import query sequences
                seqan2::StringSet<TQuerySegment, seqan2::Dependent<>> queries;
                seqan2::StringSet<seqan2::CharString> queryIDs;
                
                stellarThreadTime.input_queries_time.measure_time([&]()
                {
                    get_cart_queries(records, queries, queryIDs, thread_meta.text_out, thread_meta.text_out);
                });

                dream_stellar::_writeMoreCalculatedParams(threadOptions, threadOptions.referenceLength, queries, thread_meta.text_out);

                auto swift_index_time = stellarThreadTime.swift_index_construction_time.now();
                dream_stellar::StellarIndex<TAlphabet> stellarIndex{queries, threadOptions};
                dream_stellar::StellarSwiftPattern<TAlphabet> swiftPattern = stellarIndex.createSwiftPattern();

                // Construct index of the queries
                thread_meta.text_out << "Constructing index..." << '\n';
                stellarIndex.construct();
                thread_meta.text_out << std::endl;
                stellarThreadTime.swift_index_construction_time.manual_timing(swift_index_time);

                std::vector<size_t> disabledQueryIDs{};

                dream_stellar::StellarOutputStatistics outputStatistics{};
                bool threadFoundMatches{false};
                if (threadOptions.forward)
                {
                    auto databaseSegment = dream_stellar::_getDREAMDatabaseSegment<TAlphabet, TDatabaseSegment>
                                            (databases[threadOptions.binSequences[0]], threadOptions);
                    stellarThreadTime.forward_strand_stellar_time.measure_time([&]()
                    {
                        size_t const databaseRecordID = databaseIDMap.recordID(databaseSegment);
                        seqan2::CharString const & databaseID = databaseIDMap.databaseID(databaseRecordID);
                        // container for eps-matches
                        seqan2::StringSet<dream_stellar::QueryMatches<dream_stellar::StellarMatch<seqan2::String<TAlphabet> const,
                                                                                    seqan2::CharString> > > forwardMatches;
                        seqan2::resize(forwardMatches, length(queries));

                        constexpr bool databaseStrand = true;
                        dream_stellar::StellarComputeStatistics statistics = dream_stellar::StellarLauncher<TAlphabet>::search_and_verify
                        (
                            databaseSegment,
                            databaseID,
                            databaseStrand,
                            threadOptions,
                            swiftPattern,
                            stellarThreadTime.forward_strand_stellar_time.prefiltered_stellar_time,
                            forwardMatches
                        );

                        thread_meta.text_out << std::endl; // swift filter output is on same line
                        dream_stellar::_printDatabaseIdAndStellarKernelStatistics(threadOptions.verbose, databaseStrand, databaseID, statistics, thread_meta.text_out);

                        stellarThreadTime.forward_strand_stellar_time.post_process_eps_matches_time.measure_time([&]()
                        {
                            // forwardMatches is an in-out parameter
                            // this is the match consolidation
                            threadFoundMatches |= dream_stellar::_postproccessQueryMatches(databaseStrand, threadOptions.referenceLength, threadOptions,
                                                                    forwardMatches, disabledQueryIDs);
                        }); // measure_time

                        if (threadFoundMatches)
                        {
                            // open output files
                            std::ofstream outputFile(threadOptions.outputFile.c_str(), ::std::ios_base::out);
                            if (!outputFile.is_open())
                            {
                                std::cerr << "Could not open output file\t" << threadOptions.outputFile.c_str() << std::endl;
                                error_in_search = true;
                            }
                            stellarThreadTime.forward_strand_stellar_time.output_eps_matches_time.measure_time([&]()
                            {
                                // output forwardMatches on positive database strand
                                dream_stellar::_writeAllQueryMatchesToFile(forwardMatches, queryIDs, databaseStrand, refLen, "gff", outputFile);
                            }); // measure_time
                        }

                        outputStatistics = dream_stellar::_computeOutputStatistics(forwardMatches);
                    }); // measure_time
                }


                if (reverse)
                {
                    TDatabaseSegment databaseSegment{};
                    stellarThreadTime.reverse_complement_database_time.measure_time([&]()
                    {
                        databaseSegment = _getDREAMDatabaseSegment<TAlphabet, TDatabaseSegment>
                                            (reverseDatabases[threadOptions.binSequences[0]], threadOptions, reverse);
                    }); // measure_time

                    stellarThreadTime.reverse_strand_stellar_time.measure_time([&]()
                    {
                        size_t const databaseRecordID = reverseDatabaseIDMap.recordID(databaseSegment);
                        seqan2::CharString const & databaseID = reverseDatabaseIDMap.databaseID(databaseRecordID);
                        // container for eps-matches
                        seqan2::StringSet<dream_stellar::QueryMatches<dream_stellar::StellarMatch<seqan2::String<TAlphabet> const,
                                                                                    seqan2::CharString> > > reverseMatches;
                        seqan2::resize(reverseMatches, length(queries));

                        constexpr bool databaseStrand = false;
                        dream_stellar::StellarComputeStatistics statistics = dream_stellar::StellarLauncher<TAlphabet>::search_and_verify
                        (
                            databaseSegment,
                            databaseID,
                            databaseStrand,
                            threadOptions,
                            swiftPattern,
                            stellarThreadTime.reverse_strand_stellar_time.prefiltered_stellar_time,
                            reverseMatches
                        );

                        thread_meta.text_out << std::endl; // swift filter output is on same line
                        dream_stellar::_printDatabaseIdAndStellarKernelStatistics(threadOptions.verbose, databaseStrand, databaseID,
                                                                            statistics, thread_meta.text_out);

                        stellarThreadTime.reverse_strand_stellar_time.post_process_eps_matches_time.measure_time([&]()
                        {
                            // reverseMatches is an in-out parameter
                            // this is the match consolidation
                            threadFoundMatches |= dream_stellar::_postproccessQueryMatches(databaseStrand, threadOptions.referenceLength, threadOptions,
                                                                    reverseMatches, disabledQueryIDs);
                        }); // measure_time

                        if (threadFoundMatches)
                        {
                            // open output files
                            std::ofstream outputFile(threadOptions.outputFile.c_str(), ::std::ios_base::app);
                            if (!outputFile.is_open())
                            {
                                std::cerr << "Could not open output file\t" << threadOptions.outputFile.c_str() << std::endl;
                                error_in_search = true;
                            }
                            stellarThreadTime.reverse_strand_stellar_time.output_eps_matches_time.measure_time([&]()
                            {
                                // output reverseMatches on negative database strand
                                dream_stellar::_writeAllQueryMatchesToFile(reverseMatches, queryIDs, databaseStrand, refLen, "gff", outputFile);
                            }); // measure_time
                        }
                        outputStatistics.mergeIn(dream_stellar::_computeOutputStatistics(reverseMatches));
                    }); // measure_time
                }

                if (threadFoundMatches)
                    thread_meta.output_files.push_back(cart_queries_path.string() + ".gff");

                // Writes disabled query sequences to disabledFile.
                if (disabledQueriesFile.is_open())
                {
                    stellarThreadTime.output_disabled_queries_time.measure_time([&]()
                    {
                        // write disabled query file
                        dream_stellar::_writeDisabledQueriesToFastaFile(disabledQueryIDs, queryIDs, queries, disabledQueriesFile);
                    }); // measure_time
                }

                dream_stellar::_writeOutputStatistics(outputStatistics, threadOptions.verbose, disabledQueriesFile.is_open(), thread_meta.text_out);

                thread_meta.time_statistics.emplace_back(stellarThreadTime.milliseconds() / 1000);
                if (arguments.write_time)
                {
                    stellarThreadTime.manual_timing(current_time);
                    dream_stellar::_print_stellar_app_time(stellarThreadTime, thread_meta.text_out);
                }
                    
            }
        });
    }

    auto start = std::chrono::high_resolution_clock::now();
    // producer threads are created here
    if constexpr (stellar_only)
    {
        iterate_all_queries<TSequence>(ref_meta.seg_count, arguments, queue);
    }
    else
    {
        using index_t = decltype(index);
        raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};
        if constexpr (is_split)
        {
            iterate_split_queries<index_t, TSequence>(arguments, index, thresholder, queue, query_meta.value());
        }
        else
        {
            iterate_short_queries<index_t, TSequence>(arguments, index, thresholder, queue);
        }
    }

    queue.finish(); // Flush carts that are not empty yet
    consumerThreads.clear();
    auto end = std::chrono::high_resolution_clock::now();
    time_statistics.search_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    // merge output files and metadata from threads
    bool error_in_merge = merge_processes(arguments, time_statistics, exec_meta, var_pack);

    return error_in_search && error_in_merge;
}

}  // namespace valik::app
