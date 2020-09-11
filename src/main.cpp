#include "SG.h"
#include <string.h>
#include <unistd.h>

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1a-dirty"
#endif

int main_extendAnchor(int argc, char* argv[]);
int main_getSG(int argc, char* argv[]);
int main_mergeSG(int argc, char* argv[]);

int usage()
{
  std::cout << "\n";
  std::cout << "Program: DRAGoM (guided reference assembly of ncRNA)\n";
  std::cout << "Version: " << PACKAGE_VERSION << "\n";
  std::cout << "Contact: Ben Liu <ben_0522@ku.edu>\n\n";
  std::cout << "Usage:   DRAGoM <command> [options]\n\n";
  std::cout << "Command: \n";
  std::cout << "         getSG              generate string graph from asqg file\n";
  std::cout << "         mergeSG            connecting terminal nodes in string graph using some contigs\n";
  std::cout << "         extendAnchor       get all path from SGA string graph using anchor extension algorithm\n";
  std::cout << "\n";
  std::cout <<
  "Note: To use DRAGoM, you need to first generate assembly graph using SGA or SPAdes\n";
	return 1;
}

int main(int argc, char* argv[])
{
  clock_t t;
  t = clock();
  int ret;
  std::cout << "Welcome to DRAGoM\n";
  if (argc < 2) return usage();
  if (strcmp(argv[1], "getSG") == 0)                  ret = main_getSG(argc-1, argv+1);
  else if(strcmp(argv[1], "mergeSG") == 0)            ret = main_mergeSG(argc-1, argv+1);
  else if(strcmp(argv[1], "extendAnchor") == 0)       ret = main_extendAnchor(argc-1, argv+1);
  else
  {
    std::cout << "[main] unrecognized command "<< argv[1] << "\n";
    return usage();
  }

  return ret;
}
