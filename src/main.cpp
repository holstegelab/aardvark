#include <iostream>
#include <unistd.h>
#include <string>
#include "commands.hpp"

void printHelp(){
   printf("\nUsage:\n aardvark [command]\n");
   printf("      kc        Cyclical-shifting-aware kmer counting.\n");
   printf("      mk        Mask sequences using max-mask criteria.\n");
 }

/**
 * Main method. Parse input commands
 */
int main(int argc, char **argv){
  if(argc == 1) printHelp();
  else{
    if(std::string(argv[1]) == "kc") command_kc(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "mk") command_mk(argc - 1, &argv[1]);
    else printHelp();
  }

  return 0;
}
