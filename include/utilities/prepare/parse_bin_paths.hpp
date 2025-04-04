#pragma once

#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>

namespace valik
{

inline auto parse_bin_paths(build_arguments const & arguments, std::string const & extension = "minimiser")
{
    std::vector<std::string> minimiser_files{};
    if (arguments.bin_path.size() > 1)
    {
        for (std::filesystem::path bin_file : arguments.bin_path)
        {
            std::string bin_header = arguments.out_dir / bin_file.stem().replace_extension(extension);
            minimiser_files.push_back(bin_header);
        }
    }
    else
    {
        for (size_t bin{0}; bin < arguments.bins; bin++)
        {
            std::filesystem::path const ref_file{arguments.bin_path[0]};
            std::filesystem::path file = arguments.out_dir;
            file /= ref_file.stem();
            file += ".";
            file += std::to_string(bin);
            file += "." + extension;
            minimiser_files.emplace_back(file.string());
        }
    }

    return minimiser_files;
}

} // namespace valik
