#pragma once

#include <filesystem>
#include <unordered_map>
#include <ranges>

#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>

#include <utilities/consolidate/io.hpp>
#include <utilities/consolidate/stellar_match.hpp>

namespace valik
{

/**
 * @brief Function that removes duplicates from split Stellar search results and writes the final output file.
 *
 * @param arguments Command line arguments.
 * @return false If could not clean up intermediate results.
 */
void consolidate_matches(search_arguments const & arguments);

} // namespace valik
