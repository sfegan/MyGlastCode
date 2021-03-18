//-*-mode:c++; mode:font-lock;-*-

#include <string>
#include <iostream>
#include <fstream>

#include "Magic7Dispatcher.hpp"

int main(int argc, char** argv)
{
  std::string progname(*argv);
  argv++,argc--;

  unsigned n(0);
  Magic7PacketTool m7;
  while(argc)
    {
      n = m7.load(*argv);
      std::cerr << *argv << ": " << n << std::endl;
      argv++, argc--;
    }

  std::cerr << "Sorting... " << std::flush;
  m7.sort();
  std::cerr << std::endl;
  std::cerr << "Pruning... " << std::flush;
  n = m7.deleteDuplicates();
  std::cerr << n << std::endl;
  std::cerr << "Writing... " << std::flush;
  n = m7.save(std::cout);
  std::cerr << n << std::endl;
}
