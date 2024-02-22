#include <valik/search/search.hpp>

namespace valik::app
{

template <typename func_t>
void runtime_to_compile_time(func_t const & func, bool b1)
{
    if (b1) 
        func.template operator()<true>();
    else
        func.template operator()<false>();
}

template <typename func_t, std::same_as<bool> ...bs_t>
void runtime_to_compile_time(func_t const & func, bool b1, bs_t... bs)
{
    runtime_to_compile_time([&]<bool...compile_bs>()
    {
        if (b1)
            func.template operator()<true, compile_bs...>();
        else
            func.template operator()<false, compile_bs...>();
    }, bs...);
}

/**
 * @brief Function that loads the index and launches local or distributed search.
 *
 * @tparam compressed Interleaved Bloom Filter layout type.
 * @param arguments Command line arguments.
 */
void valik_search(search_arguments & arguments)
{
    search_time_statistics time_statistics{};

    bool failed;
    if (arguments.distribute)
    {
        runtime_to_compile_time([&]<bool is_compressed>()
        {
            failed = search_distributed<is_compressed>(arguments, time_statistics);
        }, arguments.compressed);
    }
    // Shared memory execution
    else
    {
        runtime_to_compile_time([&]<bool is_compressed, bool is_split>()
        {
            failed = search_local<is_compressed, is_split>(arguments, time_statistics);
        }, arguments.compressed, arguments.split_query);
    }

    // Consolidate matches (not necessary when searching a metagenomic database)
    auto start = std::chrono::high_resolution_clock::now();
    if (!arguments.ref_meta_path.empty())
    {
        consolidate_matches(arguments);
        const bool error_in_delete = !std::filesystem::remove(arguments.all_matches);
        if (error_in_delete)
            std::cerr << "Could not clean up intermediate file: \t" << std::string(arguments.all_matches) << '\n';
        failed = failed || error_in_delete;
    }
    auto end = std::chrono::high_resolution_clock::now();
    time_statistics.consolidation_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    if (arguments.write_time)
        write_time_statistics(time_statistics, arguments.out_file.string() + ".time", arguments);

    if (failed)
        throw std::runtime_error("valik_search failed. Run didn't complete correctly.");
}

} // namespace valik::app
