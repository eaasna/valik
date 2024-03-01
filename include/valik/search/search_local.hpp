#pragma once

#include "future"

#include <valik/search/cart_query_io.hpp>
#include <valik/search/iterate_queries.hpp>
#include <valik/search/load_index.hpp>
#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>
#include <utilities/cart_queue.hpp>
#include <utilities/consolidate/merge_processes.hpp>
#include <utilities/threshold/search_kmer_profile.hpp>

#include <stellar/database_id_map.hpp>
#include <stellar/diagnostics/print.tpp>
#include <stellar/stellar_launcher.hpp>
#include <stellar/stellar_output.hpp>
#include <stellar/io/import_sequence.hpp>
#include <stellar/utils/stellar_app_runtime.hpp>
#include <stellar3.shared.hpp>

namespace valik::app
{

/**
 * @brief Function that calls Valik prefiltering and launches parallel threads of Stellar search.
 *
 * @tparam compressed IBF layout type.
 * @tparam is_split Split query sequences.
 * @param arguments Command line arguments.
 * @param time_statistics Run-time statistics.
 * @return false if search failed.
 */
template <bool compressed, bool is_split>
bool search_local(search_arguments & arguments, search_time_statistics & time_statistics)
{
    using index_structure_t = std::conditional_t<compressed, index_structure::ibf_compressed, index_structure::ibf>;
    auto index = valik_index<index_structure_t>{};

    {
        auto start = std::chrono::high_resolution_clock::now();
        load_index(index, arguments.index_file);
        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    }

    metadata ref_meta = metadata(arguments.ref_meta_path);
    env_var_pack var_pack{};
    std::optional<metadata> query_meta;
    if (arguments.split_query)
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
            //!TODO: search profile is processed twice
            // 1. extract parameters pattern_size, max_segment_len, (<- needed before query split) and threshold
            // 2. access FNR after metadata
            std::filesystem::path search_profile_file{arguments.ref_meta_path};
            search_profile_file.replace_extension("arg");
            search_kmer_profile search_profile{search_profile_file};
            search_error_profile error_thresh = search_profile.get_error_profile(arguments.errors);     

            search_pattern pattern(arguments.errors, arguments.pattern_size);
            param_space space;
            param_set params(arguments.shape_size, arguments.threshold, space);
            filtering_request request(pattern, ref_meta, query_meta.value());
            if ((request.fpr(params) > 0.2) && (arguments.search_type != STELLAR))
                std::cerr << "WARNING: Prefiltering will be inefficient for a high error rate.\n";

            if (arguments.verbose)
            {
                std::cout.precision(3);

                if (arguments.search_type == STELLAR)
                {
                    std::cerr << "WARNING: Prefiltering will be inefficient for a high error rate. " << 
                                 "Launching exact Stellar search without prefiltering.";
                }
                else
                {
                    std::cout << "\n-----------Search parameters-----------\n";
                    std::cout << "kmer size " << std::to_string(arguments.shape_size) << '\n';
                    switch (arguments.search_type)
                    {
                        case HEURISTIC: std::cout << "heuristic "; break;
                        case LEMMA: std::cout << "k-mer lemma "; break;
                    }
                    std::cout << "threshold ";
                    std::cout << std::to_string(arguments.threshold) << '\n';

                    std::cout << "FNR " << arguments.fnr << '\n';
                    std::cout << "FPR " << request.fpr(params) << '\n';
                }
            }
        }
    }

    using TAlphabet = seqan2::Dna;
    using TSequence = seqan2::String<TAlphabet>;
    // the queue hands records over from the producer threads (valik prefiltering) to the consumer threads (stellar search) 
    auto queue = cart_queue<shared_query_record<TSequence>>{index.ibf().bin_count(), arguments.cart_max_capacity, arguments.max_queued_carts};

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

    auto bin_paths = index.bin_path();
    if (bin_paths.size() > 1 || bin_paths[0].size() > 1)
        throw std::runtime_error("Multiple reference files can not be searched in shared memory mode. "
                                 "Add --distribute argument to launch multiple distributed instances of DREAM-Stellar search.");

    stellar::stellar_runtime input_databases_time{};
    bool const databasesSuccess = input_databases_time.measure_time([&]()
    {
        std::cout << "Launching stellar search on a shared memory machine...\n";
        //!TODO: allow metagenome database
        return stellar::_importAllSequences(bin_paths[0][0].c_str(), "database", databases, databaseIDs, refLen, std::cout, std::cerr);
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
    stellar::DatabaseIDMap<TAlphabet> databaseIDMap{databases, databaseIDs};
    stellar::DatabaseIDMap<TAlphabet> reverseDatabaseIDMap{reverseDatabases, databaseIDs};

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

                thread_meta.output_files.push_back(cart_queries_path.string() + ".gff");

                stellar::StellarOptions threadOptions{};
                stellar::stellar_app_runtime stellarThreadTime{};
                threadOptions.alphabet = "dna";            // Possible values: dna, rna, protein, char
                threadOptions.verbose = true;
                threadOptions.queryFile = cart_queries_path.string();
                threadOptions.prefilteredSearch = true;
                threadOptions.referenceLength = refLen;
                threadOptions.searchSegment = true;
                auto seg = ref_meta.segment_from_bin(bin_id);
                //!TODO: deal with metagenome database
                threadOptions.binSequences.emplace_back(seg.seq_vec[0]);
                threadOptions.segmentBegin = seg.start;
                threadOptions.segmentEnd = seg.start + seg.len;
                threadOptions.numEpsilon = arguments.error_rate;
                threadOptions.epsilon = stellar::utils::fraction::from_double(threadOptions.numEpsilon).limit_denominator();
                threadOptions.minLength = arguments.pattern_size;
                threadOptions.outputFile = cart_queries_path.string() + ".gff";
                
                {
                    threadOptions.maxRepeatPeriod = arguments.maxRepeatPeriod;
                    threadOptions.minRepeatLength = arguments.minRepeatLength;
                    threadOptions.strVerificationMethod = arguments.strVerificationMethod;
                    threadOptions.xDrop = arguments.xDrop;
                    threadOptions.qgramAbundanceCut = arguments.qgramAbundanceCut;
                }

                using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
                using TQuerySegment = seqan2::Segment<seqan2::String<TAlphabet> const, seqan2::InfixSegment>;

                stellar::_writeFileNames(threadOptions, thread_meta.text_out);
                stellar::_writeSpecifiedParams(threadOptions, thread_meta.text_out);
                stellar::_writeCalculatedParams(threadOptions, thread_meta.text_out);   // calculate qGram
                thread_meta.text_out << std::endl;

                // import query sequences
                seqan2::StringSet<TQuerySegment, seqan2::Dependent<>> queries;
                seqan2::StringSet<seqan2::CharString> queryIDs;
                get_cart_queries(records, queries, queryIDs, thread_meta.text_out, thread_meta.text_out);

                stellar::_writeMoreCalculatedParams(threadOptions, threadOptions.referenceLength, queries, thread_meta.text_out);

                auto current_time = stellarThreadTime.swift_index_construction_time.now();
                stellar::StellarIndex<TAlphabet> stellarIndex{queries, threadOptions};
                stellar::StellarSwiftPattern<TAlphabet> swiftPattern = stellarIndex.createSwiftPattern();

                // Construct index of the queries
                thread_meta.text_out << "Constructing index..." << '\n';
                stellarIndex.construct();
                thread_meta.text_out << std::endl;
                stellarThreadTime.swift_index_construction_time.manual_timing(current_time);

                std::vector<size_t> disabledQueryIDs{};

                stellar::StellarOutputStatistics outputStatistics{};
                if (threadOptions.forward)
                {
                    auto databaseSegment = stellar::_getDREAMDatabaseSegment<TAlphabet, TDatabaseSegment>
                                            (databases[threadOptions.binSequences[0]], threadOptions);
                    stellarThreadTime.forward_strand_stellar_time.measure_time([&]()
                    {
                        size_t const databaseRecordID = databaseIDMap.recordID(databaseSegment);
                        seqan2::CharString const & databaseID = databaseIDMap.databaseID(databaseRecordID);
                        // container for eps-matches
                        seqan2::StringSet<stellar::QueryMatches<stellar::StellarMatch<seqan2::String<TAlphabet> const,
                                                                                    seqan2::CharString> > > forwardMatches;
                        seqan2::resize(forwardMatches, length(queries));

                        constexpr bool databaseStrand = true;
                        auto queryIDMap = stellar::QueryIDMap<TAlphabet>(queries);

                        stellar::StellarComputeStatistics statistics = stellar::StellarLauncher<TAlphabet>::search_and_verify
                        (
                            databaseSegment,
                            databaseID,
                            queryIDMap,
                            databaseStrand,
                            threadOptions,
                            swiftPattern,
                            stellarThreadTime.forward_strand_stellar_time.prefiltered_stellar_time,
                            forwardMatches
                        );

                        thread_meta.text_out << std::endl; // swift filter output is on same line
                        stellar::_printDatabaseIdAndStellarKernelStatistics(threadOptions.verbose, databaseStrand, databaseID, statistics, thread_meta.text_out);

                        stellarThreadTime.forward_strand_stellar_time.post_process_eps_matches_time.measure_time([&]()
                        {
                            // forwardMatches is an in-out parameter
                            // this is the match consolidation
                            stellar::_postproccessQueryMatches(databaseStrand, threadOptions.referenceLength, threadOptions,
                                                                    forwardMatches, disabledQueryIDs);
                        }); // measure_time

                        // open output files
                        std::ofstream outputFile(threadOptions.outputFile.c_str(), ::std::ios_base::out);
                        if (!outputFile.is_open())
                        {
                            std::cerr << "Could not open output file." << std::endl;
                            error_in_search = true;
                        }
                        stellarThreadTime.forward_strand_stellar_time.output_eps_matches_time.measure_time([&]()
                        {
                            // output forwardMatches on positive database strand
                            stellar::_writeAllQueryMatchesToFile(forwardMatches, queryIDs, databaseStrand, "gff", outputFile);
                        }); // measure_time

                    outputStatistics = stellar::_computeOutputStatistics(forwardMatches);
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
                        seqan2::StringSet<stellar::QueryMatches<stellar::StellarMatch<seqan2::String<TAlphabet> const,
                                                                                    seqan2::CharString> > > reverseMatches;
                        seqan2::resize(reverseMatches, length(queries));

                        constexpr bool databaseStrand = false;

                        stellar::QueryIDMap<TAlphabet> queryIDMap{queries};

                        stellar::StellarComputeStatistics statistics = stellar::StellarLauncher<TAlphabet>::search_and_verify
                        (
                            databaseSegment,
                            databaseID,
                            queryIDMap,
                            databaseStrand,
                            threadOptions,
                            swiftPattern,
                            stellarThreadTime.reverse_strand_stellar_time.prefiltered_stellar_time,
                            reverseMatches
                        );

                        thread_meta.text_out << std::endl; // swift filter output is on same line
                        stellar::_printDatabaseIdAndStellarKernelStatistics(threadOptions.verbose, databaseStrand, databaseID,
                                                                            statistics, thread_meta.text_out);

                        stellarThreadTime.reverse_strand_stellar_time.post_process_eps_matches_time.measure_time([&]()
                        {
                            // reverseMatches is an in-out parameter
                            // this is the match consolidation
                            stellar::_postproccessQueryMatches(databaseStrand, threadOptions.referenceLength, threadOptions,
                                                                    reverseMatches, disabledQueryIDs);
                        }); // measure_time

                        // open output files
                        std::ofstream outputFile(threadOptions.outputFile.c_str(), ::std::ios_base::app);
                        if (!outputFile.is_open())
                        {
                            std::cerr << "Could not open output file." << std::endl;
                            error_in_search = true;
                        }
                        stellarThreadTime.reverse_strand_stellar_time.output_eps_matches_time.measure_time([&]()
                        {
                            // output reverseMatches on negative database strand
                            stellar::_writeAllQueryMatchesToFile(reverseMatches, queryIDs, databaseStrand, "gff", outputFile);
                        }); // measure_time

                    outputStatistics.mergeIn(stellar::_computeOutputStatistics(reverseMatches));
                    }); // measure_time
                }

                // Writes disabled query sequences to disabledFile.
                if (disabledQueriesFile.is_open())
                {
                    stellarThreadTime.output_disabled_queries_time.measure_time([&]()
                    {
                        // write disabled query file
                        stellar::_writeDisabledQueriesToFastaFile(disabledQueryIDs, queryIDs, queries, disabledQueriesFile);
                    }); // measure_time
                }

                stellar::_writeOutputStatistics(outputStatistics, threadOptions.verbose, disabledQueriesFile.is_open(), thread_meta.text_out);

                thread_meta.time_statistics.emplace_back(stellarThreadTime.milliseconds() / 1000);
                if (arguments.write_time)
                    stellar::_print_stellar_app_time(stellarThreadTime, thread_meta.text_out);
            }
        });
    }

    auto start = std::chrono::high_resolution_clock::now();
    //!TODO: do not make thresholder if search_kind == STELLAR 
    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};
    // producer threads are created here
    if constexpr (is_split)
    {
        iterate_split_queries(arguments, index.ibf(), thresholder, queue, query_meta.value());
    }
    else
    {
        iterate_short_queries(arguments, index.ibf(), thresholder, queue);
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
