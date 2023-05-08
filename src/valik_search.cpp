#include <future>
#include <cstdlib>
#include <algorithm>

#include <valik/search/load_index.hpp>
#include <valik/search/query_record.hpp>
#include <valik/search/search_time_statistics.hpp>
#include <valik/search/prefilter_queries_parallel.hpp>
#include <valik/split/reference_metadata.hpp>
#include <valik/split/reference_segments.hpp>
#include <utilities/external_process.hpp>

#include <raptor/threshold/threshold.hpp>

#include <seqan3/io/sequence_file/output.hpp>


namespace valik::app
{

/* Creates a temporary folder in the temporary path of the OS
 *
 * \param name: a name with 'XXXXXX' at the end, e.g.: valik/call_XXXXXX
 * \return returns the name with the 'XXXXXX' replaced and the directory created
 *
 * throws if any errors occures
 */
static std::filesystem::path create_temporary_path(std::filesystem::path name) {
    if (!name.is_relative()) {
        throw std::runtime_error("Must be given a relative file");
    }
    auto path = std::filesystem::temp_directory_path() / name;
    auto path_str = path.native();
    create_directories(path.parent_path());
    auto str = std::vector<char>(path_str.size()+1, '\0'); // Must include an extra character to include a 0
    std::copy_n(path_str.data(), path_str.size(), str.data());
    auto ptr = mkdtemp(str.data());
    if (!ptr) {
        throw std::runtime_error("Could not create temporary folder: " + path_str);
    }
    return str.data();
}



//-----------------------------
//
// Setup IBF and launch multithreaded search.
//
//-----------------------------
template <bool compressed>
bool run_program(search_arguments const &arguments, search_time_statistics & time_statistics)
{
    using index_structure_t = std::conditional_t<compressed, index_structure::ibf_compressed, index_structure::ibf>;
    auto index = valik_index<index_structure_t>{};

    bool error_triggered = false; // indicates if an error happen inside this function

    {
        auto start = std::chrono::high_resolution_clock::now();
        load_index(index, arguments.index_file);
        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    }

    using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq>;
    using types = seqan3::type_list<std::string, std::vector<seqan3::dna4>>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;

    seqan3::sequence_file_input<dna4_traits, fields> fin{arguments.query_file};
    std::vector<query_record> query_records{};

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};

    // the location of bin-query fasta files can be overwritten with an environment variable
    // the $VALIK_TMP directory has to exist and write permission must be granted
    std::filesystem::path tmp_path;
    if (auto ptr = std::getenv("VALIK_TMP"); ptr != nullptr)
        tmp_path = std::string(ptr);
    else
        tmp_path = create_temporary_path("valik/stellar_call_XXXXXX");

    std::string stellar_exec = "stellar";
    if (auto ptr = std::getenv("VALIK_STELLAR"); ptr != nullptr)
        stellar_exec = std::string(ptr);

    std::string merge_exec = "cat";
    if (auto ptr = std::getenv("VALIK_MERGE"); ptr != nullptr)
        merge_exec = std::string(ptr);


    sync_out synced_out{arguments.out_file};
    auto queue = cart_queue<query_record>{index.ibf().bin_count(), arguments.cart_max_capacity, arguments.max_queued_carts};

    std::optional<reference_segments> segments;
    if (!arguments.seg_path.empty())
        segments = reference_segments(arguments.seg_path);

    std::optional<reference_metadata> ref_meta;
    if (!arguments.ref_meta_path.empty())
        ref_meta = reference_metadata(arguments.ref_meta_path, false);

    //!WORKAROUND: Stellar does not allow smaller error rates
    double er_rate = std::max((double) arguments.errors / (double) arguments.pattern_size, 0.00001);

    std::ofstream text_out(arguments.out_file.string() + ".out");


    std::mutex mutex;
    std::unordered_map<size_t, size_t> bin_count;
    std::vector<std::string> output_files;

    auto consumerThreads = std::vector<std::jthread>{};
    for (size_t i = 0; i < arguments.threads; ++i)
    {
        consumerThreads.emplace_back(
        [&]() {
            // this will block until producer threads have added carts to queue
            for (auto next = queue.dequeue(); next; next = queue.dequeue())
            {
                auto & [bin_id, records] = *next;
		if (records.empty())
			continue;
                std::unique_lock g(mutex);
                std::filesystem::path path = tmp_path / std::string("query_" + std::to_string(bin_id) + "_" + std::to_string(bin_count[bin_id]++) + ".fasta");

                output_files.push_back(path.string() + ".gff");
                g.unlock();

                {
                    seqan3::sequence_file_output fout{path, fields{}};

                    for (auto & record : records)
                    {
                        sequence_record_type sequence_record{std::move(record.sequence_id), std::move(record.sequence)};
                        fout.push_back(sequence_record);
                    }
                }

                std::vector<std::string> process_args{};
                if (arguments.write_time)
                {
                    std::filesystem::path time_path =  path.string() + std::string(".gff.time");
                    process_args.insert(process_args.end(), {"/usr/bin/time", "-o", std::string(time_path), "-f", "\"%e\t%M\t%x\t%C\""});
                }
                process_args.push_back(stellar_exec);

                if (segments && ref_meta)
                {
                    auto ref_len = ref_meta->total_len;
                    auto seg = segments->segment_from_bin(bin_id);
                    process_args.insert(process_args.end(), {index.bin_path()[0][0], std::string(path),
                                                            "--referenceLength", std::to_string(ref_len),
                                                            "--sequenceOfInterest", std::to_string(seg.ref_ind),
                                                            "--segmentBegin", std::to_string(seg.start),
                                                            "--segmentEnd", std::to_string(seg.start + seg.len)});
                }
                else
                {
                    process_args.insert(process_args.end(), {index.bin_path()[bin_id][0], std::string(path)});
                }

                process_args.insert(process_args.end(), {"-e", std::to_string(er_rate),
                                                        "-l", std::to_string(arguments.pattern_size),
                                                        "-o", std::string(path) + ".gff"});

                auto start = std::chrono::high_resolution_clock::now();
                external_process process(process_args);
                auto end = std::chrono::high_resolution_clock::now();
                time_statistics.cart_processing_times.push_back(0.0 + std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count());

                text_out << process.cout();
                text_out << process.cerr();

                if (process.status() != 0) {
                    std::unique_lock g(mutex); // make sure that our output is synchronized
                    std::cerr << "error running VALIK_STELLAR\n";
                    std::cerr << "call:";
                    for (auto args : process_args) {
                        std::cerr << " " << args;
                    }
                    std::cerr << '\n';
                    std::cerr << process.cerr() << '\n';
                    error_triggered = true;
                }
            }

        });
    }

    for (auto &&chunked_records : fin | seqan3::views::chunk((1ULL << 20) * 10))
    {
        query_records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        for (auto && query_record: chunked_records)
            query_records.emplace_back(std::move(query_record.id()), std::move(query_record.sequence()));

        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        start = std::chrono::high_resolution_clock::now();
        prefilter_queries_parallel(index.ibf(), arguments, query_records, thresholder, queue);
        end = std::chrono::high_resolution_clock::now();
        time_statistics.prefilter_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    }
    queue.finish(); // Flush carts that are not empty yet
    consumerThreads.clear();

    std::vector<std::string> merge_process_args{merge_exec};
    for (auto & path : output_files)
        merge_process_args.push_back(path);
    external_process merge(merge_process_args);

    auto check_external_process_success = [&](std::vector<std::string> const & proc_args, external_process const & proc)
    {
        if (proc.status() != 0) {
            std::cerr << "External process failed\n";
            std::cerr << "call:";
            for (auto args : proc_args) {
                std::cerr << " " << args;
            }
            std::cerr << '\n';
            std::cerr << proc.cerr() << '\n';
            return true;
        }
        else
            return error_triggered;
    };

    if (arguments.write_time)
    {
        std::vector<std::string> merge_time_files{"cat"};
        for (auto & path : output_files)
            merge_time_files.push_back(path + ".time");
        external_process merge_time(merge_time_files);
        error_triggered = check_external_process_success(merge_time_files, merge_time);

        std::ofstream time_out(arguments.out_file.string() + std::string(".time"));
        time_out << merge_time.cout();
    }

    error_triggered = check_external_process_success(merge_process_args, merge);

    std::ofstream matches_out(arguments.out_file);
    matches_out << merge.cout();
    return error_triggered;
}

void valik_search(search_arguments const & arguments)
{
    search_time_statistics time_statistics{};

    bool failed;
    if (arguments.compressed)
        failed = run_program<seqan3::data_layout::compressed>(arguments, time_statistics);
    else
        failed = run_program<seqan3::data_layout::uncompressed>(arguments, time_statistics);

    if (arguments.write_time)
        write_time_statistics(time_statistics, arguments);

    if (failed) {
        throw std::runtime_error("valik_search failed. Run didn't complete correctly.");
    }
}

} // namespace valik::app
