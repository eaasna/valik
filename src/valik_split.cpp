#include <valik/argument_parsing/shared.hpp>
#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>
#include <valik/split/write_seg_sequences.hpp>

namespace valik::app
{

//-----------------------------
//
// Divide reference or query database into partially overlapping segments.
//
//-----------------------------
void valik_split(split_arguments const & arguments)
{
    metadata meta(arguments);
    meta.to_file(arguments.meta_out);

    if (arguments.write_ref)
        write_reference_segments(meta, arguments.meta_out);
    if (arguments.write_query)
        write_query_segments(meta, arguments.meta_out);

    metadata meta_deserialised(arguments.meta_out);
}

} // namespace valik::app
