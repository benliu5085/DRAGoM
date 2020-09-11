#include <unistd.h>
#include "SG.h"

std::string CPU;
std::string DIR_CDHIT;

int extendAnchor_usage()
{
  std::cout << "\n";
  std::cout << "Usage:   DRAGoM extendAnchor [options] <config> <string_graph.fq> <anchor1> <anchor2> ... \n\n";
  std::cout << "options: \n";
  std::cout << "         -m INT     maximum length when building a path, [100]\n";
  std::cout << "                    m < 0 to turn off hard cutoff, extend the path to match model length\n";
  std::cout << "         -d BOOL    anchor masking parameter, off by default\n";
  std::cout << "                    speed up at the possible cost of losing sensitivity\n";
  std::cout << "         -c INT     start compress when <INT> paths generated\n";
  std::cout << "         -t INT     use <INT> CPUs to compress paths\n";
  std::cout << "         -f STR     directory to cd-hit-est\n";
  std::cout << "\n";
	return 1;
}

int main_extendAnchor(int argc, char* argv[])
{

  if(argc <= 1) return extendAnchor_usage();
  int max_length = 100, c, batchsize = 2000000;
  bool aggressive = false;

  while ((c = getopt(argc, argv, "dc:m:t:f:")) >= 0)
  {
    if (c == 'm')  max_length = std::stoi(optarg);
    else if (c == 'd')  aggressive = true;
    else if (c == 'c')  batchsize = std::stoi(optarg);
    else if (c == 't')  CPU = optarg;
    else if (c == 'f')  DIR_CDHIT = optarg;
    else                return extendAnchor_usage();
  }

  std::string fq_name = argv[optind];
  optind++;
  StrGraph* p_g = new StrGraph();
  p_g->readFqFile(fq_name);
  std::string cm_name, path_name;
  while(optind < argc)
  {
    cm_name = argv[optind];
    path_name = cm_name.substr(0, cm_name.rfind('.')) + ".path.fa";
    std::cout << "Now extending " << cm_name << "\n";
    p_g->DirectedDFS(cm_name, max_length, path_name, aggressive, batchsize);
    optind++;
  }

  delete p_g;
  return 0;
}
