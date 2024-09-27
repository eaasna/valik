#pragma once

#include <memory>
#include <valik/split/metadata.hpp>

namespace valik
{

template <typename sequence_t>
struct shared_reference_record
{
    std::shared_ptr<sequence_t> underlying_sequence;
    metadata::segment_stats segment;

    shared_reference_record(std::shared_ptr<sequence_t> const & seq_ptr, metadata::segment_stats && seg) : 
                            underlying_sequence(seq_ptr), segment(seg) { }
};

} // namespace valik
