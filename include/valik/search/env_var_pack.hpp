#pragma once

#include <cstdlib>
#include <unistd.h>
#include <filesystem>
#include <vector>
#include <algorithm>

#include <utilities/external_process.hpp>
#include <valik/search/execution_metadata.hpp>
#include <valik/search/search_time_statistics.hpp>
#include <valik/shared.hpp>

namespace valik
{

struct env_var_pack
{
    std::filesystem::path tmp_path;
    std::string stellar_exec{"stellar"};

    env_var_pack()
    {
        // the location of bin-query fasta files can be overwritten with an environment variable
        if (auto ptr = std::getenv("VALIK_TMP"); ptr != nullptr)
        {
            tmp_path = std::string(ptr);
            std::filesystem::file_status s = status(tmp_path);
            std::filesystem::perms p = s.permissions();
            bool is_writable = (std::filesystem::perms::none != (p & std::filesystem::perms::owner_write)) |
                               (std::filesystem::perms::none != (p & std::filesystem::perms::group_write)) |
                               (std::filesystem::perms::none != (p & std::filesystem::perms::others_write));

            if (!exists(tmp_path) | !is_directory(s) | !is_writable )
                throw std::runtime_error("Directory $VALIK_TMP=" + std::string(ptr) + " must exist and write permission must be granted");
        }

        else
            tmp_path = create_temporary_path("valik/stellar_call_XXXXXX");

        if (auto ptr = std::getenv("VALIK_STELLAR"); ptr != nullptr)
            stellar_exec = std::string(ptr);
    }

    /* Creates a temporary folder in the temporary path of the OS
    *
    * \param name: a name with 'XXXXXX' at the end, e.g.: valik/call_XXXXXX
    * \return returns the name with the 'XXXXXX' replaced and the directory created
    *
    * throws if any errors occurs
    */
    static std::filesystem::path create_temporary_path(std::filesystem::path name)
    {
        if (!name.is_relative())
        {
            throw std::runtime_error("Must be given a relative file");
        }
        auto path = std::filesystem::temp_directory_path() / name;
        auto path_str = path.native();
        create_directories(path.parent_path());
        auto str = std::vector<char>(path_str.size()+1, '\0'); // Must include an extra character to include a 0
        std::copy_n(path_str.data(), path_str.size(), str.data());
        auto ptr = mkdtemp(str.data());
        if (!ptr)
        {
            throw std::runtime_error("Could not create temporary folder: " + path_str);
        }
        return str.data();
    }
};

}
