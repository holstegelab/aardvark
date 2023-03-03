#include "commands.hpp"
#include "kc.hpp"
#include "cxxopts.hpp"
#include <vector>
#include <string>

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_kc_parser(int argc, char** argv){
  cxxopts::Options options(argv[0]);
  try{
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <FASTA[.GZ]>...")
      .add_options()
      ("k, kmer-length", "Kmer-length", cxxopts::value<int>())
      ("m, max-k", "Iteratively count kmers up to this value if > 0", cxxopts::value<int>()->default_value("0"));
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.emplace_back(i);
    if(inputs.empty()) std::cout << options.help();
    else {
      int k = result["kmer-length"].as<int>();
      if(k < 1 || k > 32){
        std::cerr << "ERROR: kmer-length must be > 0 and < 32" << std::endl;
        exit(1);
      }
      int m = result["max-k"].as<int>();
      if(m == 0) m = k + 1;
      else if(m < 0 || m <= k) {
        std::cerr << "ERROR: max-k must be > 0 and > k" << std::endl;
        exit(1);
      }
      kc(inputs, (uint32_t)k, (uint32_t)m);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " <<  e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_kc(int argc, char **argv){
  //parse CLI arguments
  command_kc_parser(argc, argv);
  return 0;
}
