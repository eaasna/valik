#pragma once

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include <iostream>
#include <fstream>

namespace sliding_window
{

class reference_metadata
{
    public:
        struct sequence_stats 
        {
            std::string id;             
            size_t len;        
        };

        size_t total_len;
        std::vector<sequence_stats> sequences;

        reference_metadata(std::string filename) 
        { 
            using traits_type = seqan3::sequence_file_input_default_traits_dna;
            seqan3::sequence_file_input<traits_type> fin{filename};

            total_len = 0;
            for (auto & record : fin)
            {
                sequence_stats seq;
                seq.id = record.id();
                seq.len = record.sequence().size();
                sequences.push_back(seq);
                total_len += seq.len;
            }
        }

        void to_file(std::string filepath)
        {
            std::ofstream out_file;
            out_file.open(filepath);
            out_file << "ID" << '\t' << "LEN" << '\n';
            for (sequence_stats & seq : sequences)
            {
                out_file << seq.id << '\t' << seq.len << '\n';
            }

            out_file << "TOTAL LEN: " << total_len << '\n';
            out_file.close();
        }

};

} // namespace sliding_window