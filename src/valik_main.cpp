#include <seqan3/argument_parser/all.hpp>

#include <valik/argument_parsing/split.hpp>
#include <valik/argument_parsing/build.hpp>
#include <valik/argument_parsing/search.hpp>
#include <valik/argument_parsing/shared.hpp>
#include <valik/argument_parsing/top_level.hpp>
#include <valik/valik.hpp>

int main(int argc, char ** argv)
{
    try
    {
        sharg::parser top_level_parser{"valik", argc, argv, sharg::update_notifications::off, {"build", "search"}};
        valik::app::init_top_level_parser(top_level_parser);

        valik::app::try_parsing(top_level_parser);

        sharg::parser & sub_parser = top_level_parser.get_sub_parser();
        if (sub_parser.info.app_name == std::string_view{"valik-build"})
            valik::app::run_build(sub_parser);
        if (sub_parser.info.app_name == std::string_view{"valik-search"})
            valik::app::run_search(sub_parser);
    }
    catch(std::exception const& e)
    {
        std::cerr << "[Error] " << e.what() << '\n';
        std::exit(-1);
    }
    catch(...)
    {
        std::cerr << "[Error] unknown exception type\n";
        std::exit(-1);
    }

    return 0;
}
