#pragma once

#include <valik/argument_parsing/shared.hpp>
#include <valik/build/build.hpp>
#include <valik/split/metadata.hpp>
#include <utilities/prepare/compute_bin_size.hpp>

namespace valik::app
{

void init_build_parser(sharg::parser & parser, build_arguments & arguments);
void run_build(sharg::parser & parser);

} // namespace valik::app
