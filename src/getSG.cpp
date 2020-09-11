#include "SG.h"
bool renameASQG(char* fname)
{
  std::ifstream fin;
  fin.open(fname);
  std::string newname = fname;
  std::string outname1 = newname.substr(0,newname.find_last_of('.')) + ".rename.asqg";
  std::string outname2 = newname.substr(0,newname.find_last_of('.')) + ".namemap.txt";
  std::ofstream fout1(outname1.c_str());
  std::ofstream fout2(outname2.c_str());
  if (!fin.is_open())
  {
    std::cout << "Error: unable to open " << newname << "\n";
    return false;
  }
  else
  {
    std::string line;
    StringVector    p_header;
    String2Integer  p_readID;

    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    while(std::getline(fin, line))
    {
      if(line[0] == 'V')
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();            // cont[0], VT
        fout1 << *it << "\t";
        it++;                           // cont[1], read name
        int id = p_readID.size();
        p_readID[*it] = id;
        fout1 << id << "\t";
        fout2 << *it << "\n";
        it++;                           // cont[2], sequence
        fout1 << *it << "\n";
      }
      else if(line[0] == 'E')
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();            // cont[0], ED
        fout1 << *it << "\t";
        it++;                           // cont[1], read1 name
        int source_id = p_readID[*it];
        fout1 << source_id << " ";
        it++;                           // cont[2], read2 name
        int target_id = p_readID[*it];
        fout1 << target_id;
        it++;
        while(it != tokens.end())
        {
          fout1 << " " << *it;
          it++;
        }
        fout1 << "\n";
      }
      else
      {
        fout1 << line << "\n";
      }
    }
  }
  fin.close();
  fout1.close();
  fout2.close();
  return true;
}

int main_getSG(int argc, char* argv[])
{
  if(argc != 2)
  {
    std::cout << "\n";
    std::cout << "Usage:   DRAGoM getSG <asqg> \n\n";
    return 1;
  }
  else
  {
    std::string newname = argv[1];
    std::string outname1 = newname.substr(0,newname.find_last_of('.')) + ".rename.asqg";
    if(renameASQG(argv[1]))
    {
      StrGraph* p_g = new StrGraph();
      p_g->readAsqgFile(outname1);
      p_g->CondenseGraph();
      p_g->writeGraph(outname1);
      delete p_g;
      return 0;
    }
    else return 1;
  }
}
