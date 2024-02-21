#include <utilities/threshold/io.hpp>

namespace valik
{

std::filesystem::path fn_filename()
{
    auto space = param_space();
    std::string outfile = "fn_err_conf_e" + std::to_string(space.max_errors) + "_t" + std::to_string(space.max_thresh) + 
                                        "_l" + std::to_string(space.max_len) + "_k" + std::to_string(std::get<0>(space.kmer_range)) + "_" + 
                                        std::to_string(std::get<1>(space.kmer_range)) + ".tsv";
    return outfile;
}


/**
 * @brief Precalculate and write out FN configuration count tables.
*/
void precalc_fn_confs(std::vector<kmer_loss> & attr_vec)
{
    auto space = param_space();
    auto outfile = fn_filename();
    for (auto k = std::get<0>(space.kmer_range); k <= std::get<1>(space.kmer_range); k++)
    {
        kmer_loss attr(k);
        attr_vec.push_back(attr);
        std::ofstream outstr;
        outstr.open(outfile, std::ios_base::app);
        attr.serialize(outstr);
    }
}

/**
 * @brief Deserialize kmer_loss for a single k-mer size.
*/
kmer_loss deserialize_kmer_loss(std::string const & kmer_attr_str)
{
    std::istringstream matrix_str(kmer_attr_str);
    std::string line;
    std::getline(matrix_str, line);
    uint8_t k{0};
    try
    {
        k = std::stoi(line.substr(line.find("k=") + 2));
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << "\nCan't convert " << line.substr(line.find("k=") + 2)  << "to k-mer size\n";
    }
    
    std::vector<std::vector<std::vector<uint64_t>>> matrix;
    std::vector<std::vector<uint64_t>> table;
    for (; std::getline(matrix_str, line); )
    {
        bool is_first_thresh_row{false};
        try
        {
            is_first_thresh_row = (stoi(line.substr(line.find("t=") + 2)) == 1);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << "\nCan't convert " << line.substr(line.find("t=") + 2) << " to threshold\n";
        }
        
        std::istringstream row_str(line);
        if (line.find("t=") == std::string::npos)
        {
            std::string cell;
            std::vector<uint64_t> row;
            for (; std::getline(row_str, cell, '\t'); )
            {
                try
                {
                    row.push_back(std::stoull(cell));
                }
                catch(const std::exception& e)
                {
                    std::cerr << e.what() << "\nCan't convert cell value\t" << cell << '\n';
                }
            }
            table.push_back(row);
        }
        else if (is_first_thresh_row)
        {
            continue;
        }
        else
        {
            matrix.push_back(table);
            table.clear();
        }
    }
    
    matrix.push_back(table);
    return kmer_loss(k, matrix);
}

/**
 * @brief Read precalculated FN error configuration count tables across all k-mer sizes. 
*/
bool read_fn_confs(std::vector<kmer_loss> & attr_vec)
{
    auto infile = fn_filename();
    if (!std::filesystem::exists(infile))
        return false;

    std::ifstream instr(infile);
    std::string line;
    std::string kmer_attr_str;
    bool past_first(false);
    while (std::getline(instr, line)) 
    {
        if (line.find("k=") != std::string::npos)
        {
            if (past_first)
            {
                attr_vec.push_back(deserialize_kmer_loss(kmer_attr_str));
                kmer_attr_str.clear();
            }
            past_first = true;
        }
        kmer_attr_str += line + '\n';
    }

    attr_vec.push_back(deserialize_kmer_loss(kmer_attr_str));
    return true;
}

}   // namespace valik
