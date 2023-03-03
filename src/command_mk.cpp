#include "commands.hpp"
#include "mk.hpp"
#include "cxxopts.hpp"
#include <vector>
#include <string>

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_mk_parser(int argc, char** argv){
  cxxopts::Options options(argv[0]);
  try{
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <FASTA[.GZ]>...")
      .add_options()
      ("m, mask-seqs", "Sequences to mask", cxxopts::value<std::string>());
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.emplace_back(i);
    if(inputs.empty()) std::cout << options.help();
    else {
      const std::string s = result["mask-seqs"].as<std::string>();
      mk(inputs, s);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_mk(int argc, char **argv){
  //parse CLI arguments
  command_mk_parser(argc, argv);
  return 0;
}
