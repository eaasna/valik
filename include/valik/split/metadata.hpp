#pragma once

#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include <utilities/threshold/param_set.hpp>
#include <utilities/threshold/search_pattern.hpp>
#include <utilities/threshold/basics.hpp>
#include <valik/shared.hpp>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <ranges>
#include <sstream>

#include <cereal/archives/binary.hpp> 
#include <cereal/types/vector.hpp>

namespace valik
{

/**
 * @brief Function that removes leading and trailing whitespace.
*/
template <typename id_t>
void trim_fasta_id(id_t & id)
{
    std::string whitespace = " \t\n\r\v";
    auto first_valid = id.find_first_not_of(whitespace);
    if (first_valid == std::string::npos)
        throw std::runtime_error{"Sequence name can not be empty."};

    std::istringstream iss(id);
    iss >> id;
}

/**
 * @brief Struct that stores the metadata for a split database.
 *  \param total_len    Total database length.
 *  \param seq_count    Number of sequences.
 *  \param seg_count    Database was divided into this many segments.
 *  \param sequences    Collection of database sequences.
 *  \param ibf_fpr  FPR of a k-mer query in the IBF.
 *  \param default_seg_len  Default length of a segment that is dynamically updated.
 *  \param segments     Collection of database segments.
 */
struct metadata
{
    /** !\brief a metadata struct that represents a sequence file.
     *
     * \param id    Numerical file id.    
     * \param path  Input fasta file path.
     */
    struct sequence_file
    {
        size_t id;
        std::string path;

        sequence_file() noexcept = default;
        sequence_file(sequence_file const &) noexcept = default;
        sequence_file & operator=(sequence_file const &) noexcept = default;
        sequence_file & operator=(sequence_file &&) noexcept = default;
        ~sequence_file() noexcept = default;
        
        sequence_file(size_t const i, std::string const p) : id(i), path(p) { }

        template <class Archive>
        void serialize(Archive & archive)
        {
            archive(id, path);
        }
    };

    /** !\brief a metadata struct that represents a single sequence.
     *
     * \param file_id   Numerical file id (sequence_file::id).
     * \param id        The FASTA id.
     * \param ind       0-based index in the input FASTA file.
     * \param len       Sequence length.
     */
    struct sequence_stats
    {
        size_t file_id;
        std::string id;
        size_t ind;
        uint64_t len;

        sequence_stats() noexcept = default;
        sequence_stats(sequence_stats const &) noexcept = default;
        sequence_stats & operator=(sequence_stats const &) noexcept = default;
        sequence_stats & operator=(sequence_stats &&) noexcept = default;
        ~sequence_stats() noexcept = default;
        
        sequence_stats(size_t const seq_file_id, std::string const fasta_id, size_t const fasta_ind, uint64_t const seq_length) :
                       file_id(seq_file_id), id(fasta_id), ind(fasta_ind), len(seq_length) {}

        template <class Archive>
        void serialize(Archive & archive)
        {
            archive(file_id, id, ind, len);
        }
    };

    struct length_order
    {
        inline bool operator() (sequence_stats const & left, sequence_stats const & right)
        {
            return (left.len < right.len);
        }
    };

    /** !\brief a struct that represents a single segment of a reference or query database.
     *
     * All indices and positions are 0-based.
     *
     *  \param id           Numerical segment id.
     *  \param seq_vec      List of sequences (numerical ids corresponding to sequence_stats::ind) associated with this segment.
     *  \param start        Segment start position in sequence if segment consists of a single subsequence. 0 for metagenome bin.
     *  \param len          Segment length.
     */
    struct segment_stats
    {
        size_t id;
        std::vector<size_t> seq_vec{};
        uint64_t start{0};
        uint64_t len;

        segment_stats() noexcept = default;
        segment_stats(segment_stats const &) noexcept = default;
        segment_stats & operator=(segment_stats const &) noexcept = default;
        segment_stats & operator=(segment_stats &&) noexcept = default;
        ~segment_stats() noexcept = default;

        segment_stats(size_t const ind, uint64_t const s, uint64_t const l) : start(s), len(l)
        {
            seq_vec.push_back(ind);
        }

        segment_stats(size_t const & i, std::vector<size_t> & ind_vec, uint64_t const l) : id(i), len(l) 
        { 
            seq_vec = std::move(ind_vec);
        }

        std::string unique_id()
        {
            std::string str_id{};
            for (auto seq_ind : seq_vec)
            {
                str_id += std::to_string(seq_ind);
                str_id += "_";
            }
            str_id += std::to_string(start);
            str_id += "_";
            str_id += std::to_string(len);
            return str_id;
        }

        template <class Archive>
        void serialize(Archive & archive)
        {
            archive(id, seq_vec, start, len);
        }
    };

    struct fasta_order
    {
        inline bool operator() (sequence_stats const & left, sequence_stats const & right)
        {
            return (left.ind < right.ind);
        }

        inline bool operator() (segment_stats const & left, segment_stats const & right)
        {
            if (left.seq_vec.size() > 1 || right.seq_vec.size() > 1)
                throw std::runtime_error("Can't order sets of sets of sequences.");
            return (left.seq_vec[0] < right.seq_vec[0]);
        }
    };

    uint64_t total_len{0};
    size_t seq_count;
    size_t seg_count;
    size_t pattern_size;
    float ibf_fpr;

    std::vector<sequence_file> files;
    std::vector<sequence_stats> sequences;
    std::vector<segment_stats> segments;

    private:    
        size_t default_seg_len;

        /**
         * @brief Function that scans over a sequence file to extract metadata.
         *
         * @param db_path Path to input file.
         */
        void scan_database_file(std::string const & db_file)
        {
            using traits_type = seqan3::sequence_file_input_default_traits_dna;
            files.emplace_back(0, db_file);
            seqan3::sequence_file_input<traits_type> fin{db_file};   // single input file
            size_t fasta_ind = sequences.size();
            for (auto & record : fin)
            {
                trim_fasta_id(record.id());
                sequence_stats seq(0, record.id(), fasta_ind, record.sequence().size());    // there is only a single sequence file
                total_len += seq.len;
                sequences.push_back(seq);
                fasta_ind++;
            }
            std::stable_sort(sequences.begin(), sequences.end(), length_order());
        }

        void add_segment(size_t const seq_id, uint64_t const s, uint64_t const l)
        {
            segments.emplace_back(seq_id, s, l);
        }

        void add_segment(size_t const i, std::vector<size_t> & bin_seq_ids, uint64_t const l)
        {
            segments.emplace_back(i, bin_seq_ids, l);
        }

        /**
         * @brief Function that scans over a metagenome database to extract sequences and segments.
         *
         * @param bin_path Database paths.
         */
        void scan_metagenome_bins(std::vector<std::vector<std::string>> const & bin_path)
        {
            using traits_type = seqan3::sequence_file_input_default_traits_dna;
            size_t file_id{0};
            for (std::vector<std::string> bin_files : bin_path)
            {
                uint64_t bin_len{0};
                std::vector<size_t> bin_seq_ids;
                for (std::string bin_file : bin_files)
                {
                    files.emplace_back(file_id, bin_file);
                    seqan3::sequence_file_input<traits_type> fin{bin_file};
                    size_t fasta_ind = sequences.size();
                    for (auto & record : fin)
                    {
                        trim_fasta_id(record.id());
                        sequence_stats seq(file_id, record.id(), fasta_ind, record.sequence().size());
                        total_len += seq.len;
                        bin_len += seq.len;
                        bin_seq_ids.push_back(fasta_ind);
                        sequences.push_back(seq);
                        fasta_ind++;
                    }
                    file_id++;
                }
                add_segment(segments.size(), bin_seq_ids, bin_len); 
            }
        }

        /**
         * @brief Function that splits the database into n partially overlapping segments.
         *
         * @param n Number of segments.
         * @param overlap Length of overlap between adjacent segments.
         * @param seq_it Iterator to first sequence of sufficient length.
         */
        template <typename it_t>
        void make_exactly_n_segments(size_t const & n, size_t const & overlap, it_t & seq_it)
        {
            size_t remaining_db_len = total_len;
            for (auto it = seq_it; it != sequences.end(); it++)
            {
                auto seq = *it;
                size_t start = 0;
                if (seq.len <= default_seg_len * 1.5)
                {
                    // database sequence is single segment
                    add_segment(seq.ind, start, seq.len);
                    remaining_db_len -= seq.len;
                }
                else
                {
                    // sequences that are contained in a single segment might not have the exact segment length
                    // dynamically update segment length to divide the rest of the remaining database as equally as possible among the chosen number of segments
                    size_t remaining_seg_count = n - segments.size();
                    size_t updated_seg_len = std::round((float) remaining_db_len / remaining_seg_count);
                    size_t segments_per_seq = std::round( (double) seq.len / (double) updated_seg_len);

                    if (segments_per_seq == 1)
                        add_segment(seq.ind, start, seq.len);
                    else
                    {
                        size_t actual_seg_len = std::ceil(((float) seq.len - overlap) / segments_per_seq);

                        // divide database sequence into multiple segments
                        add_segment(seq.ind, 0, actual_seg_len + overlap);

                        for (start += actual_seg_len; start + actual_seg_len + overlap < seq.len - overlap; start += actual_seg_len)
                        {
                            add_segment(seq.ind, start, actual_seg_len + overlap);
                        }
                        add_segment(seq.ind, start, seq.len - start);
                    }
                    remaining_db_len -= seq.len;
                }
            }

            if (segments.size() != n)
            {
                throw std::runtime_error("Database was split into " + std::to_string(segments.size()) +
                                         " instead of " + std::to_string(n) + " segments.");
            }
        }

        /**
         * @brief Function that splits the database into partially overlapping segments of roughly equal length.
         *
         * @param overlap Length of overlap between adjacent segments.
         * @param seq_it Iterator to first sequence of sufficient length.
         */
        template <typename it_t>
        void make_equal_length_segments(size_t const & overlap, it_t & seq_it)
        {
            for (auto it = seq_it; it != sequences.end(); it++)
            {
                auto seq = *it;
                size_t start = 0;
                if (seq.len <= default_seg_len * 1.5)
                {
                    // database sequence is single segment
                    add_segment(seq.ind, start, seq.len);
                }
                else
                {
                    size_t segments_per_seq = std::round( (double) seq.len / (double) default_seg_len);

                    if (segments_per_seq == 1)
                        add_segment(seq.ind, start, seq.len);
                    else
                    {
                        size_t actual_seg_len = std::ceil(((float) seq.len - overlap) / segments_per_seq);

                        // divide database sequence into multiple segments
                        add_segment(seq.ind, 0, actual_seg_len + overlap);

                        for (start += actual_seg_len; start + actual_seg_len + overlap < seq.len - overlap; start += actual_seg_len)
                        {
                            add_segment(seq.ind, start, actual_seg_len + overlap);
                        }
                        add_segment(seq.ind, start, seq.len - start);
                    }
                }
            }
        }

        /**
         * @brief Function that scans over the collection of sequences and greedily assigns segment positions.
         *
         * @param seg_count_in Chosen number of segments.
         * @param n Actual number of segments.
         * @param overlap Length of overlap between adjacent segments.
         */
        template <typename arg_t>
        void scan_database_sequences(arg_t const & arguments)
        {
            default_seg_len = total_len / arguments.seg_count + 1;
            if (default_seg_len <= arguments.pattern_size)
            {
                throw std::runtime_error("Segments of length " + std::to_string(default_seg_len) + "bp can not overlap by " +
                                         std::to_string(arguments.pattern_size) + "bp.\nDecrease the overlap or the number of segments.");
            }

            size_t len_lower_bound = default_seg_len / 10;
            // Check how many sequences are discarded for being too short
            //!TODO: gather short sequences
            auto first_long_seq = std::find_if(sequences.begin(), sequences.end(), [&](auto const & seq){return (seq.len > len_lower_bound);});
            size_t discarded_short_sequences = first_long_seq - sequences.begin();

            for (auto it = sequences.begin(); it < first_long_seq; it++)
            {
                seqan3::debug_stream << "Sequence: " << (*it).id << " is too short and will be skipped.\n";
                total_len -= (*it).len;
            }

            if (arguments.seg_count < (sequences.size() - discarded_short_sequences))
            {
                throw std::runtime_error("Can not split " + std::to_string(sequences.size()) +
                                         " sequences into " + std::to_string(arguments.seg_count) + " segments.");
            }

            if constexpr (std::is_same<arg_t, split_arguments>::value)
                make_exactly_n_segments(arguments.seg_count, arguments.pattern_size, first_long_seq);
            else
                make_equal_length_segments(arguments.pattern_size, first_long_seq);

            std::stable_sort(sequences.begin(), sequences.end(), fasta_order());
            std::stable_sort(segments.begin(), segments.end(), fasta_order());
            for (size_t i = 0; i < segments.size(); i++)
                segments[i].id = i;
        }

    public:

        template <typename arg_t>
        void update_segments_for_distributed_stellar(arg_t arguments)
        {
            segments.clear();
            default_seg_len = total_len / arguments.seg_count + 1;
            if (default_seg_len <= arguments.pattern_size)
            {
                throw std::runtime_error("Segments of length " + std::to_string(default_seg_len) + "bp can not overlap by " +
                                         std::to_string(arguments.pattern_size) + "bp.\nDecrease the overlap or the number of segments.");
            }

            auto seq_begin = sequences.begin();
            make_equal_length_segments(arguments.pattern_size, seq_begin);
            seg_count = segments.size();

            std::stable_sort(sequences.begin(), sequences.end(), fasta_order());
            std::stable_sort(segments.begin(), segments.end(), fasta_order());
            for (size_t i = 0; i < segments.size(); i++)
                segments[i].id = i;
                
        }

        /**
         * @brief Constructor that scans a sequence database to create a metadata struct.
        */
        metadata(split_arguments const & arguments)
        {
            ibf_fpr = arguments.fpr;
            if (arguments.metagenome)
            {
                scan_metagenome_bins(arguments.bin_path);
            }
            else
            {
                scan_database_file(arguments.bin_path[0][0]);
                scan_database_sequences(arguments);
            }

            seq_count = sequences.size();
            seg_count = segments.size();
            pattern_size = arguments.pattern_size;
        }

        metadata(search_arguments & arguments)
        {
            scan_database_file(arguments.query_file);
            if (!arguments.manual_parameters)
            {
                //!TODO: what is a suitable value for this magic const
                //!TODO: switch between prefiltered and stellar_only search?
                if (total_len > (arguments.max_segment_len * 10))
                    arguments.seg_count = std::round(total_len / (arguments.max_segment_len - arguments.pattern_size));
                else
                    arguments.seg_count = std::max(arguments.seg_count, (uint32_t) sequences.size() * 2);
            }
            else
            {
                arguments.seg_count = arguments.seg_count_in;
            }

            scan_database_sequences(arguments);
            if (arguments.manual_parameters && (segments.size() != arguments.seg_count))
                seqan3::debug_stream << "[Warning] Database was split into " << segments.size() << " instead of " << arguments.seg_count << " segments.\n";

            seq_count = sequences.size();
            seg_count = segments.size();
            pattern_size = arguments.pattern_size;
        }

        /**
         * @brief Constructor that deserializes a metadata struct from file.
         */
        metadata(std::filesystem::path const & filepath)
        {
            load(filepath);
        }

        /**
         * @brief Function that returns the numerical index of a sequence based on its fasta ID.
         *
         * @param string_id Fasta ID.
         */
        inline size_t ind_from_id(std::string const & string_id) const
        {
            auto it = std::find_if(sequences.begin(), sequences.end(), [&](const sequence_stats& seq) { return seq.id == string_id;});
            if (it == sequences.end())
                throw seqan3::validation_error{"Sequence metadata does not contain sequence " + string_id + " from Stellar output."};
            else
                return (*it).ind;
        }

        /**
         * @brief Function that returns the segment corresponding to a numerical ID.
         *
         * @param id
         * @return segment
         */
        segment_stats segment_from_bin(size_t const id) const
        {
            if (segments.size() <= id)
                throw std::runtime_error{"Segment " + std::to_string(id) + " index out of range."};

            return segments[id];
        }

        /**
         * @brief Function that returns the slice of segments corresponding to a sequence.
         *
         * @param ind Index of sequence.
         */
        auto segments_from_ind(size_t const & ind)
        {
            if (sequences.size() <= ind)
                throw std::runtime_error{"Sequence " + std::to_string(ind) + " index out of range."};

            return segments | std::views::filter([ind](segment_stats const & seg) 
            {
                return (std::find(seg.seq_vec.begin(), seg.seq_vec.end(), ind) != seg.seq_vec.end());
            });
        }

        /**
         * @brief Serialize the metadata struct.
         *
         * @param filepath Output file path.
         */
        void save(std::filesystem::path const & filepath) const
        {
            std::ofstream os(filepath, std::ios::binary);
            cereal::BinaryOutputArchive archive(os);
            archive(total_len, pattern_size, files, sequences, segments, ibf_fpr);
        }
      
        /**
         * @brief Deserialise the metadata struct.
         *
         * @param filepath Input file path.
         */
        void load(std::filesystem::path const & filepath)
        {
            std::ifstream is(filepath, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(total_len, pattern_size, files, sequences, segments, ibf_fpr);
            seq_count = sequences.size();
            seg_count = segments.size();
        }

        std::string to_string()
        {
            std::stringstream out_str;
            for (sequence_stats const & seq : sequences)
                out_str << seq.id << '\t' << seq.ind << '\t' << seq.len << '\n';

            out_str << "$\n";

            for (size_t seg_id{0}; seg_id < segments.size(); seg_id++)
            {
                segment_stats seg = segments[seg_id];
                out_str << seg_id << '\t';
                for (size_t ind : seg.seq_vec) 
                    out_str << ind << '\t';
                out_str << seg.start << '\t' << seg.len << '\n';
            }

            out_str << "$\n";

            return out_str.str();
        }

        double segment_length_stdev()
        {
            auto segment_lengths = segments | std::views::transform([](segment_stats const & seg){return seg.len;});

            double sum = std::accumulate(segment_lengths.begin(), segment_lengths.end(), 0.0);
            double mean = sum / segment_lengths.size();
            double sq_sum = std::inner_product(segment_lengths.begin(), segment_lengths.end(), segment_lengths.begin(), 0.0);
            double stdev = std::sqrt(sq_sum / segments.size() - mean * mean);

            return stdev;
        }

        // coefficient of variation
        double segment_length_cv()
        {
            auto segment_lengths = segments | std::views::transform([](segment_stats const & seg){return seg.len;});

            double sum = std::accumulate(segment_lengths.begin(), segment_lengths.end(), 0.0);
            double mean = sum / segment_lengths.size();
            double sq_sum = std::inner_product(segment_lengths.begin(), segment_lengths.end(), segment_lengths.begin(), 0.0);
            double stdev = std::sqrt(sq_sum / segments.size() - mean * mean);

            return stdev / mean;
        }

        /**
        * @brief Probability of a k-mer appearing spuriously in a bin.
        */
        double kmer_spurious_match_prob(uint8_t const kmer_size) const
        {
            return std::min(1.0, expected_kmer_occurrences(std::round(total_len / (double) seg_count), kmer_size));
        }

        /**
        * @brief The probability of at least threshold k-mers matching spuriously between a query pattern and a reference bin.
        */
        double pattern_spurious_match_prob(param_set const & params) const
        {
            double fpr{1};
            double p = ibf_fpr + kmer_spurious_match_prob(params.k);
            const double precision{1e-9};
            /*
            For parameters 
        
            s reference size
            b number of bins
            k kmer size
            l min match length  (equal to segment overlap)
            t threshold 

            find the probability of t or more k-mers matching spuriously between a pattern of length l and a reference bin.   
        
            The probability of a pattern exceeding the threshold by chance is equal to 1 - P(0 k-mers match) - 
                                                                                       P(1 k-mer matches) - 
                                                                                       ... - 
                                                                                       P(t - 1 k-mers match).

            Given p which is the probability of a single k-mer appearing spuriously in the reference bin and n = (l - k + 1) k-mers per pattern 
            the probability of 0 k-mers matching spuriously is P(0 k-mers match) = (1 - p)^n 
            and
            the probability of 1 k-mer out of n matching is P(1 k-mer matches) = (n take 1) * p^1 * (1 - p)^(n - 1).
            */

            size_t kmers_per_pattern = pattern_size - params.k + 1;
            for (uint8_t matching_kmer_count{0}; matching_kmer_count < params.t; matching_kmer_count++)
            {
                if (fpr < precision) 
                    break;
                fpr -= combinations(matching_kmer_count, kmers_per_pattern) * 
                                    pow(p, matching_kmer_count) * 
                                    pow(1 - p, kmers_per_pattern - matching_kmer_count);
            }

            return std::max(0.0, fpr - precision);
        }
};

} // namespace valik
