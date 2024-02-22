#pragma once

#include <utilities/threshold/find.hpp>
#include <utilities/threshold/io.hpp>
#include <valik/argument_parsing/shared.hpp>
#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>
#include <valik/split/write_seg_sequences.hpp>


namespace valik::app
{

void valik_split(split_arguments & arguments);

} // namespace valik::app
