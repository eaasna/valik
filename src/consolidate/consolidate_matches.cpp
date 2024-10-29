#include <utilities/consolidate/consolidate_matches.hpp>

namespace valik
{

/**
 * @brief Function that truncates the fasta id if it contains whitespace.
*/
std::string truncate_fasta_id(std::string const & id)
{
    auto first_whitespace = id.find_first_of(" ");
    if (first_whitespace == std::string::npos)
        return id;

    return id.substr(0, first_whitespace);
}

void consolidate_matches(search_arguments const & arguments)
{
    auto ref_meta = metadata(arguments.ref_meta_path);
    auto matches = read_alignment_output<stellar_match>(arguments.all_matches, ref_meta);
    //auto before_duplicate_removal = matches.size();
    std::sort(matches.begin(), matches.end(), std::greater<stellar_match>());
    matches.erase( std::unique( matches.begin(), matches.end() ), matches.end() );

    //seqan3::debug_stream << before_duplicate_removal << '\t' << matches.size() << '\n';
    // <query_ind, <refs>>
    std::unordered_multimap<std::string, std::string> overabundant_pairs{}; 
    std::unordered_set<std::string> disabled_queries{};

    // <query_id, <ref_id, match_count>>
    std::unordered_map<std::string, std::unordered_map<std::string, size_t>> total_match_counter{};
    decltype(matches) consolidated_matches{};
    
    for (auto & match : matches)
    {
        if ( total_match_counter[match.qname][match.dname] < arguments.disableThresh )
        {
            total_match_counter[match.qname][match.dname]++;
        }
    }
    
    // for <query, ref> pairs that do not appear often return all matches
    for (auto & match : matches)
    {
        size_t query_match_count{0};
        for (auto & it : total_match_counter[match.qname])
            query_match_count += it.second;

        bool is_disabled = query_match_count >= arguments.disableThresh;
        bool is_overabundant = total_match_counter[match.qname][match.dname] > arguments.numMatches;
         
        if (!is_overabundant && !is_disabled)
            consolidated_matches.emplace_back(match);
        else if (is_disabled)
            disabled_queries.emplace(match.qname);
        else
            overabundant_pairs.emplace(match.qname, match.dname);
    }

    // for <query, ref> pairs that appear often return arguments.numMatches longest matches
    for (auto & it : overabundant_pairs)
    {
        std::vector<stellar_match> overabundant_matches{};
        auto is_overabundant_match = [&](auto & m)
        {
            return ((m.qname == it.first) && (m.dname == it.second));
        };

        for (auto & m : matches | std::views::filter(is_overabundant_match))
        {
            overabundant_matches.push_back(m);
        }

        std::sort(overabundant_matches.begin(), overabundant_matches.end(), stellar_match::length_order()); // sort in order of increasing length
        consolidated_matches.insert(consolidated_matches.end(), overabundant_matches.end() - arguments.numMatches, overabundant_matches.end());    
    }

    // debug
    if (arguments.verbose && !overabundant_pairs.empty())
    {
        seqan3::debug_stream << "Overabundant pairs\n";
        for (auto & it : overabundant_pairs)
        {
            seqan3::debug_stream << it.first << '\t' << it.second << '\n'; 
        }
    }
    
    // write out fasta file of disabled queries
    if (disabled_queries.size() > 0)
    {
        using input_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>>;
        input_file_t fin{arguments.query_file};
        
        seqan3::sequence_file_output fdisabled(arguments.disabledQueriesFile);
        
        for (auto & record : fin)
        {
            if (std::find(disabled_queries.begin(), disabled_queries.end(), truncate_fasta_id(record.id())) != disabled_queries.end())
                fdisabled.push_back(record);
        }

        // debug
        if (arguments.verbose)
            seqan3::debug_stream << "Disabled " << disabled_queries.size() << " queries.\n";
    }

    write_alignment_output<stellar_match>(arguments.out_file, consolidated_matches);
}

}  // namespace valik
