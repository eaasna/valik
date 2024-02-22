#pragma once

#include <utilities/threshold/basics.hpp>
#include <utilities/threshold/filtering_request.hpp>
#include <utilities/threshold/io.hpp>

#include <valik/search/query_record.hpp>
#include <valik/search/search_distributed.hpp>
#include <valik/search/search_local.hpp>
#include <valik/search/search_time_statistics.hpp>
#include <valik/shared.hpp>

#include <utilities/consolidate/consolidate_matches.hpp>

namespace valik::app
{

void valik_search(search_arguments & arguments);

} // namespace valik::app
