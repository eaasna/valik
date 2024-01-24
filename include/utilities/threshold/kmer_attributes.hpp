#pragma once

#include <utilities/threshold/shared.hpp>

namespace valik
{

/**
 * @brief For a kmer size k consider how the FNR changes for various values of the parameters t, e and l. 
 *   
 * @param k Chosen k-mer size.
 * @param fn_conf_counts Number of configurations of e errors that do not retain at least t shared kmers after mutation.
*/
struct kmer_attributes
{
    using row_t = std::vector<uint64_t>;
    using table_t = std::vector<row_t>;
    using mat_t = std::vector<table_t>;
    
    const size_t k;
    const mat_t fn_conf_counts;  // false negative configuration counts

    /**
     * @brief For a maximum error count find the number of error configurations 
     *        that do not retain at least the threshold number of k-mers after mutation.
    */
    mat_t count_err_conf_below_thresh()
    {
        auto space = param_space();
        mat_t matrix;
        for (size_t t = 1; t <= space.max_thresh; t++)
        {
            table_t table;
            for (size_t e = 0; e <= space.max_errors; e++)
            {
                row_t row;
                for(size_t l = 0; l <= space.max_len; l++)
                {       
                    int shared_base_count = l - e;
                    if ((int) (k + t) - 1 > shared_base_count) // all error distributions destroy enough k-mers
                        row.push_back(total_err_conf_count(e, l));
                    else if (e == 0)
                        row.push_back(0);
                    else
                    {
                        row.push_back(row.back() + table.back()[l - 1]);

                        size_t band = l - k - t;
                        row.back() -= table.back()[band + t - 1];

                        for (size_t i = 1; i < t; i++)
                        {
                            row.back() -= matrix[matrix.size() - i][e - 1][band + t - i - 1];
                            row.back() += matrix[matrix.size() - i][e - 1][band + t - i];
                        }
                    }
                }
                table.push_back(row);
            }
            matrix.push_back(table);
        }
        return matrix;
    }

    kmer_attributes(size_t const kmer_size) : 
                k(kmer_size), 
                fn_conf_counts(count_err_conf_below_thresh()) { }

    kmer_attributes(size_t const kmer_size, mat_t const & matrix) : k(kmer_size), fn_conf_counts(matrix) {}

    /**
     * @brief False negative rate for a parameter set.
    */
    double fnr_for_param_set(filtering_request const & request, param_set const & params) const
    {
        return fn_conf_counts[params.t - 1][request.e][request.l] / (double) request.total_conf_count();
    }

    /**
     * @brief Score of the objective function for a parameter set. Smaller values are better.
    */
    double score(filtering_request const & request, param_set const & params) const
    {
        return fnr_for_param_set(request, params) + request.fpr(k) / (double) params.t;
    }

    /**
     * @brief Stream out flattened 3D matrix of false negative configuration counts.
    */
    template <typename str_t>
    void serialize(str_t & outstr)
    {
        outstr << "k=" << k << '\n'; 
        size_t t = 1;
        for (auto threshold_table : fn_conf_counts)
        {   
            outstr << "t=" << t << '\n';
            for (auto error_row : threshold_table)
            {
                for (auto cell : error_row)
                    outstr << cell << '\t';
                outstr << '\n';
            }
            t++;
        }
    }
};

}   // namespace valik
