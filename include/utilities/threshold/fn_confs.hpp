#pragma once

#include <utilities/threshold/kmer_loss.hpp>

namespace valik
{

/**
 * @brief Gather estimates for the number of false negatives across a range of k-mer lengths, error rates, thresholds and match lengths. 
 *   
 * @param space Parameter space.
*/
struct fn_confs
{
    param_space space;

    /**
     * @brief Serialization file path.
     */
    std::filesystem::path fn_filename()
    {
        std::string outfile = "fn_err_conf_e" + std::to_string(space.max_errors) + "_t" + std::to_string(space.max_thresh) + 
                                                "_l" + std::to_string(space.max_len) + "_k" + 
                                                std::to_string(space.min_k()) + "_" + 
                                                std::to_string(space.max_k()) + ".bin";
        return outfile;
    }

    private:
        std::vector<kmer_loss> attr_vec;

    public:
        /**
         * @brief Serialize the false negative error configurations.
         *
         * @param filepath Output file path.
         */
        void save(std::filesystem::path const & filepath) const
        {
            std::ofstream os(filepath, std::ios::binary);
            cereal::BinaryOutputArchive archive(os);
            archive(space, attr_vec);
        }

        /**
        * @brief Deserialise the false negative error configurations.
        *
        * @param filepath Input file path.
        */
        void load(std::filesystem::path const & filepath)
        {
            std::ifstream is(filepath, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(space, attr_vec);
        }

    private:
        /**
         * @brief Precalculate and write out FN configuration count tables.
        */
        void precalc_fn_confs()
        {
            auto outfile = fn_filename();
            for (auto k = space.min_k(); k <= space.max_k(); k++)
            {
                kmer_loss attr(k, space);
                attr_vec.push_back(attr);
            }
            save(outfile);
        }

        /**
         * @brief Read precalculated FN error configuration count tables across all k-mer sizes. 
        */
        bool try_read_fn_confs()
        {
            auto infile = fn_filename();
            if (!std::filesystem::exists(infile))
                return false;

            load(infile);
            return true;
        }

    public:
        fn_confs(param_space const & s) : space(s)
        {
            if (!try_read_fn_confs())
                precalc_fn_confs();
        }

        const kmer_loss & get_kmer_loss(uint8_t const k) const
        {
            if (k < space.min_k() || k > space.max_k())
            {
                throw std::runtime_error("k-mer length " + std::to_string(k) + " is out of range [" + 
                      std::to_string(space.min_k()) + ", " + std::to_string(space.max_k()) + "]");
            }
            
            return attr_vec[k - space.min_k()];
        }

        const kmer_loss & get_kmer_loss(seqan3::shape const shape)
        {
            uint8_t k = shape.size();
            if (k < space.min_k() || k > space.max_k())
            {
                throw std::runtime_error("k-mer length " + std::to_string(k) + " is out of range [" + 
                      std::to_string(space.min_k()) + ", " + std::to_string(space.max_k()) + "]");
            }

            //!TODO: adjust precalculated attribute vectors to account for gapped shapes
            // currently hack to copy the ungapped attribute without adjusting threshold 
            attr_vec[k - space.min_k()].shape = shape;
            return attr_vec[k - space.min_k()];
        }
};

} // namespace valik
