#pragma once

#include <valik/shared.hpp>
#include <valik/split/database_metadata.hpp>

#include <utilities/consolidate/io.hpp>

namespace valik
{

/**
 * @brief Function that removes duplicates from split Stellar search results and writes the final output file.
 *
 * @param arguments Command line arguments.
 */
void consolidate_matches(search_arguments const & arguments);

} // namespace valik