#include <valik/search/load_index.hpp>
#include <valik/search/query_record.hpp>
#include <valik/search/search_time_statistics.hpp>
#include <valik/search/write_output_file_parallel.hpp>

#include <raptor/threshold/threshold.hpp>

namespace valik::app
{

//-----------------------------
//
// Setup IBF and launch multithreaded search.
//
//-----------------------------
template <bool compressed>
void run_program(search_arguments const &arguments, search_time_statistics & time_statistics)
{
    using index_structure_t = std::conditional_t<compressed, index_structure::ibf_compressed, index_structure::ibf>;
    auto index = valik_index<index_structure_t>{};

    auto cereal_worker = [&]()
    {
        auto start = std::chrono::high_resolution_clock::now();
        load_index(index, arguments.index_file);
        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    };

    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{arguments.query_file};
    std::vector<query_record> query_records{};

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};

    sync_out synced_out{arguments.out_file};

    // ==========================================
    // Create sync_out for each bin-query path.
    // ==========================================
    std::vector<sync_out> bin_query_writers{};
    const std::filesystem::path bin_query_file = arguments.bin_query_file;
    std::ifstream istrm{bin_query_file};
    std::string line;

    while (std::getline(istrm, line))
    {
        if (!line.empty())
        {
            std::filesystem::path bin_query_path = line;
            sync_out bin_query_out{bin_query_path};
            bin_query_writers.emplace_back(std::move(bin_query_out));
        }
    } // sync_out can not be const; arguments are passed as const

    size_t record_id = 0u;
    for (auto &&chunked_records : fin | seqan3::views::chunk((1ULL << 20) * 10))
    {
        query_records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        for (auto && query_record: chunked_records)
        {
            query_records.emplace_back(record_id, std::move(query_record.id()), std::move(query_record.sequence()));
            ++record_id;
        }
        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        start = std::chrono::high_resolution_clock::now();

        //TODO: pass index and overlap instead of ibf and all parameters
        write_output_file_parallel(index.ibf(), arguments, query_records, thresholder, synced_out, bin_query_writers, time_statistics);
        end = std::chrono::high_resolution_clock::now();
        time_statistics.compute_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    }
}

void valik_search(search_arguments const & arguments)
{
    search_time_statistics time_statistics{};

    if (arguments.compressed)
        run_program<seqan3::data_layout::compressed>(arguments, time_statistics);
    else
        run_program<seqan3::data_layout::uncompressed>(arguments, time_statistics);

    if (arguments.write_time)
        write_time_statistics(time_statistics, arguments);

    return;
}

} // namespace valik::app
