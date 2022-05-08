#pragma once

#include <filesystem>

#include <valik/shared.hpp>

namespace valik::app
{

struct search_time_statistics
{
    double index_io_time{0.0};
    double reads_io_time{0.0};
    double bin_queries_io_time{0.0};
    double compute_time{0.0};
};

inline void write_time_statistics(search_time_statistics const & time_statistics, search_arguments const & arguments)
{
    std::filesystem::path file_path{arguments.out_file};
    file_path += ".time";
    std::ofstream file_handle{file_path};
    file_handle << "IBF I/O\tReads I/O\tBin-queries I/O\tCompute\n";
    file_handle << std::fixed
                << std::setprecision(2)
                << time_statistics.index_io_time << '\t'
                << time_statistics.reads_io_time << '\t'
                << time_statistics.bin_queries_io_time << '\t'
                << time_statistics.compute_time;

}

} // namespace valik::app
