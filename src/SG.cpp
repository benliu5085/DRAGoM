#include "SG.h"
#include "unionFind.hpp"
extern std::string CPU;
extern std::string DIR_CDHIT;

/***************************** for string graph *******************************/
bool StrGraph::readFqFile(std::string& filename)
{
  std::cout << "start loading Fq\n";
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  std::ifstream fin;
  // fin.open(config.IDIR_GRAPH.c_str());	// open file
  fin.open(filename.c_str());	// open file
  if (!fin.is_open())	    // fail to open
  {
    std::cout << "Error: " << filename << " doesnot exist!\n";
    return false;
  }
  else				// succeed to open
  {
    std::string line;
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    int cnt_E = 0;
    while(std::getline(fin, line))
    {
      if(line[0] == '>')
      {
        line.erase(0,1);
        int colon = line.find(':');
        if(colon == std::string::npos)  // no tail
        {
          std::string header = line.substr(0,colon);
          int colon = line.find(',');
          std::string path_info_str = header.substr(colon+1);
          // p_header.push_back(path_info_str);

          boost::tokenizer<boost::char_separator<char> > tokens(header, sep);
          it = tokens.begin();
          int node_id = std::stoi(*it);
          ++it;
          int read_id = -1;
          IntegerVector read_on_this_node, loc_on_this_node;
          for(int i = 0 ; it != tokens.end(); ++it, ++i)
          /* format: read_id, location, read_length */
          {
            if(i % 3 == 0)        read_on_this_node.push_back(std::stoi(*it));
            else if( i % 3 == 1)  loc_on_this_node.push_back(std::stoi(*it));
            else                  p_read_length[*(read_on_this_node.end()-1)] = std::stoi(*it);
          }
          p_read_on_node.push_back(read_on_this_node);
          p_pos_on_node.push_back(loc_on_this_node);

          BoostSTRVertex v_source;
          if(vertex.count(node_id) == 0)
          {
            STRVertexType node(node_id);
            v_source = boost::add_vertex(node, *p_graph_);
            vertex[node_id] = v_source;
          }
          else
          {
            v_source = vertex[node_id];
          }
        }
        else
        {
          std::string header = line.substr(0,colon);
          std::string tail = line.substr(colon+1);
          int colon = line.find(',');
          std::string path_info_str = header.substr(colon+1);
          // p_header.push_back(path_info_str);

          boost::tokenizer<boost::char_separator<char> > tokens(header, sep);
          it = tokens.begin();
          int node_id = std::stoi(*it);
          ++it;
          IntegerVector read_on_this_node, loc_on_this_node;
          for(int i = 0 ; it != tokens.end(); ++it, ++i)
          /* format: read_id, location, read_length */
          {
            if(i % 3 == 0)        read_on_this_node.push_back(std::stoi(*it));
            else if( i % 3 == 1)  loc_on_this_node.push_back(std::stoi(*it));
            else                  p_read_length[*(read_on_this_node.end()-1)] = std::stoi(*it);
          }
          p_read_on_node.push_back(read_on_this_node);
          p_pos_on_node.push_back(loc_on_this_node);

          BoostSTRVertex v_source, v_target;
          if(vertex.count(node_id) == 0)
          {
            STRVertexType node(node_id);
            v_source = boost::add_vertex(node, *p_graph_);
            vertex[node_id] = v_source;
          }
          else
          {
            v_source = vertex[node_id];
          }

          boost::tokenizer<boost::char_separator<char> > tokens1(tail, sep);
          for(it = tokens1.begin(); it != tokens1.end(); ++it)
          {
            int target_id = std::stoi(*it);
            if(vertex.count(target_id) > 0)
            {
              v_target = vertex[target_id];
            }
            else
            {
              STRVertexType node(target_id);
              v_target = boost::add_vertex(node, *p_graph_);
              vertex[target_id] = v_target;
            }

            std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_target, *p_graph_);
            if(e_search.second)
            {
              (*p_graph_)[e_search.first].sid_ = cnt_E;
              cnt_E++;
            }
            else
            {
              std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
            }
          }
        }
      }
      else
      {
        p_seq.push_back(line);
      }
    }
  }
  p_order = p_seq.size();
  ff_stringGraph = true;
  fin.clear();
  fin.close();
  /** timer **/
  t = clock() - t;
  std::cout << "load fastq takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  return true;
}

unsigned int goodMapping(std::string& CIAGR, int flag)
{
  if((flag&4) == 4)                 return 0x80000000;  // unmap
  unsigned int match = 0;
  std::string temp = "";
  std::string symbols = "";
  for(int i = 0; i < CIAGR.size(); i++)
  {
    if(CIAGR[i] >= '0' && CIAGR[i] <= '9')    temp += CIAGR[i];
    else
    {
        symbols += CIAGR[i];
        if(CIAGR[i] == 'M') match += std::stoi(temp);
        temp = "";
    }
  }
  if(match < 100)                   return 0x20000000;  // mapping length too short
  if(symbols.size() > 2)            return 0x20000000;  // too fragmented
  if (symbols.size() == 2)
  {
    if(((flag&16) == 16) && symbols[1] == 'M')
                                    return 0x40000000;  // tail unmap
    if(((flag&16) == 0)  && symbols[1] != 'M')
                                    return 0x40000000;  // tail unmap
    return (0x10000000+match);
  }
  return match;
}

bool goodJunction(BooleanVector& coverage, int pos_i, int pos_j)
{
  int b = pos_i, e = pos_j;
  if(pos_i > pos_j) {    b = pos_j; e = pos_i;  }
  int cnt_c = 0;
  for(int i = b; i <= e; i++)
  {
    if(!coverage[i])
      cnt_c++;
    if(cnt_c > int(0*(e-b+1))) return false;
  }
  return true;
}

bool thisIsOrphat(std::string& header)
{
  boost::char_separator<char> sep(",");
  boost::tokenizer<boost::char_separator<char> > tokens(header, sep);
  boost::tokenizer<boost::char_separator<char> >::iterator it;
  it = tokens.begin();
  ++it; // discard read_id
  int i;
  for(i = 0 ; it != tokens.end(); ++it, ++i) if(i > 3) return false;
  return (i == 3);
}

void StrGraph::mergeSG(std::string& sam_name, std::string& contig_name)
{
  /*********/
  std::ifstream fin;
  /*********/
  clock_t t;
  /*********/
  /* load contigs */
  /*********/
  t = clock();
  /*********/
  fin.open(contig_name.c_str());
  StringVector contigs;
  String2Integer contigs_id;
  if (!fin.is_open())	    // fail to open
  {
    std::cout << "Error: " << contig_name << " doesnot exist!\n";
  }
  else				// succeed to open
  {
    std::string line;
    std::string name;
    while(std::getline(fin, line))
    {
      if(line[0]=='>')
      {
        name = line.substr(1);
        contigs_id[name] = contigs_id.size()-1;
      }
      else
      {
        if(contigs.size() < contigs_id.size())  contigs.push_back("");
        contigs[contigs.size()-1] += line;
      }
    }
  }
  fin.clear();
  fin.close();
  /*********/
  t = clock() - t;
  std::cout << "load contigs takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /*********/
  /* load terminals */
  /*********/
  t = clock();
  /*********/
  std::vector<BooleanVector> taken_2D;
  for(auto x : contigs) taken_2D.push_back(BooleanVector(x.size(), false));
  std::vector<BooleanVector> coverage_2D;
  for(auto x : contigs) coverage_2D.push_back(BooleanVector(x.size(), false));
  std::vector<std::map<int, int>> f_terminals = std::vector<std::map<int, int>>(contigs.size());
  std::vector<std::map<int, int>> f_match = std::vector<std::map<int, int>>(contigs.size());
  fin.open(sam_name.c_str());	// open file
  if (!fin.is_open())	    // fail to open
  {
    std::cout << "Error: " << sam_name << " doesnot exist!\n";
  }
  else				// succeed to open
  {
    std::string line;
    boost::char_separator<char> sep("\t");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    while(std::getline(fin, line))
    {
      if(line[0] != '@')
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();  // edge name
        int node_id = std::stoi((*it));
        std::string header = (*it);
        bool ff_orphant_read = thisIsOrphat( header );
        ++it;                 // flag
        int flag = std::stoi((*it));
        ++it;                 // reference name
        int ref_id = contigs_id[(*it)];
        ++it;                 // pos
        std::string pos_str = (*it);
        ++it;                 // MAPQ
        ++it;                 // CIAGR
        std::string CIAGR = (*it);
        BoostSTRVertex node = vertex[node_id];

        // unsigned int mapping_quality = goodMapping(CIAGR, flag);
        // if(ff_orphant_read)
        // {
        //   if(mapping_quality != 0)
        //   {
        //     int read_length = p_seq[node_id].size();
        //     int pos = std::stoi(pos_str);
        //     for(int i = 0; i < read_length; ++i)
        //     {
        //       coverage_2D[ref_id][pos+i-1] = true;
        //       taken_2D[ref_id][pos+i-1] = true;
        //     }
        //   }
        // }
        // else
        // {
        //   int pos = std::stoi(pos_str);
        //   for(int i = 0; i < mapping_quality; ++i) taken_2D[ref_id][pos+i-1] = true;
        //   if((mapping_quality != 0) && ((boost::in_degree(node, *p_graph_)*boost::out_degree(node, *p_graph_)) == 0))
        //   {
        //     if((flag&16) == 0)
        //     {
        //       f_terminals[ref_id][pos] = node_id;
        //       f_match[ref_id][node_id] = mapping_quality;
        //     }
        //   }
        // }
        unsigned int mapping_quality = goodMapping(CIAGR, flag);
        if(ff_orphant_read)
        {
          if((mapping_quality&0xF0000000) == 0x00000000)
          {
            int read_length = p_seq[node_id].size();
            int pos = std::stoi(pos_str);
            for(int i = 0; i < read_length; ++i)
            {
              coverage_2D[ref_id][pos+i-1] = true;
              taken_2D[ref_id][pos+i-1] = true;
            }
          }
        }
        else
        {
          int pos = std::stoi(pos_str);
          int match_length = (mapping_quality&0x0FFFFFFF);
          for(int i = 0; i < match_length; ++i) taken_2D[ref_id][pos+i-1] = true;
          if(((mapping_quality&0xE0000000) == 0x00000000) && ((boost::in_degree(node, *p_graph_)*boost::out_degree(node, *p_graph_)) == 0))
          {
            if((flag&16) == 0)
            {
              f_terminals[ref_id][pos] = node_id;
              f_match[ref_id][node_id] = match_length;
            }
          }
        }

      }
    }
  }
  fin.clear();
  fin.close();
  /*********/
  t = clock() - t;
  std::cout << "load terminals takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /*********/

  /* building DisjointSet */
  /*********/
  t = clock();
  /*********/
  DisjointSet* ds = new DisjointSet();
  for(auto it_e = boost::edges(*p_graph_).first; it_e != boost::edges(*p_graph_).second; ++it_e)
  {
    BoostSTRVertex head = boost::source(*it_e, *p_graph_);
    BoostSTRVertex tail = boost::target(*it_e, *p_graph_);
    ds->unionNode((*p_graph_)[head].rid_, (*p_graph_)[tail].rid_);
  }
  /*********/
  t = clock() - t;
  std::cout << "construct DisjointSet takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /*********/

  /*********/
  t = clock();
  /*********/
  for(int ref_id = 0; ref_id < f_terminals.size(); ref_id++)
  {
    auto this_terminal = f_terminals[ref_id];
    auto it1 = this_terminal.begin();
    auto it2 = this_terminal.begin();
    std::string seq = contigs[ref_id];

    // take orphant seq
    int p = 0, q = 0;
    while(p < seq.size())
    {
      while((p < seq.size()) && taken_2D[ref_id][p]) p++;
      if(p - q > 100)
      {
        std::string new_seq = seq.substr(q, p-q);
        IntegerVector read_on_this_node, loc_on_this_node;
        read_on_this_node.push_back(1);
        loc_on_this_node.push_back(1);

        int new_node_id = p_seq.size();
        p_seq.push_back(new_seq);
        p_read_on_node.push_back(read_on_this_node);
        p_pos_on_node.push_back(loc_on_this_node);

        BoostSTRVertex new_node;
        if(vertex.count(new_node_id) == 0)
        {
          STRVertexType node(new_node_id);
          new_node = boost::add_vertex(node, *p_graph_);
          vertex[new_node_id] = new_node;
        }
        else
        {
          new_node = vertex[new_node_id];
        }
      }
      p++;
      q = p;
    }

    // merge
    if(it2 != this_terminal.end())
    {
      // before first
      int pos1 = it1->first;
      int node1 = it1->second;
      int s = pos1-2;
      while((s>=0) && coverage_2D[ref_id][s]) s--;
      if(s != (pos1-2))
      {
        int first_read_id = p_read_on_node[node1][0];
        int len = pos1 + p_read_length[first_read_id] - s - 2;
        std::string new_seq = seq.substr(s+1, len);
        IntegerVector read_on_this_node, loc_on_this_node;
        read_on_this_node.push_back(first_read_id);
        loc_on_this_node.push_back(pos1 - s + 1);

        int new_node_id = p_seq.size();
        p_seq.push_back(new_seq);
        p_read_on_node.push_back(read_on_this_node);
        p_pos_on_node.push_back(loc_on_this_node);

        BoostSTRVertex new_node;
        if(vertex.count(new_node_id) == 0)
        {
          STRVertexType node(new_node_id);
          new_node = boost::add_vertex(node, *p_graph_);
          vertex[new_node_id] = new_node;
        }
        else
        {
          new_node = vertex[new_node_id];
        }

        BoostSTRVertex v_target = vertex[node1];
        std::pair<BoostSTREdge, bool> e_search2 = boost::add_edge(new_node, v_target, *p_graph_);
        if(!e_search2.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
      }

      it2++;
      for(; it2 != this_terminal.end(); it1++, it2++)
      {
        int pos1 = it1->first;
        int node1 = it1->second;
        int pos2 = it2->first;
        int node2 = it2->second;
        int end = pos1 + f_match[ref_id][node1] - 2;              // coordinate of last nt of node1 on contig
        if(ds->find(node1) != ds->find(node2)) // node1 and node2 should not be connected!!!
        {
          int last_read_id = *(p_read_on_node[node1].end()-1);
          int first_read_id = p_read_on_node[node2][0];
          int l1 = p_read_length[last_read_id];
          int l2 = p_read_length[first_read_id];
          int y = end + 1 - l1;
          int z = pos2 + l2 -2;
          int len = z - y + 1;

          if(len == l1)
          {
            BoostSTRVertex v_source = vertex[node1];
            BoostSTRVertex v_target = vertex[node2];
            std::pair<BoostSTREdge, bool> e_search1 = boost::add_edge(v_source, v_target, *p_graph_);
            if(!e_search1.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
          }
          else
          {
            std::string new_seq;
            IntegerVector read_on_this_node, loc_on_this_node;
            if(len > l1)
            {
              new_seq = seq.substr(y, len);
              read_on_this_node.push_back(last_read_id);
              read_on_this_node.push_back(first_read_id);
              loc_on_this_node.push_back(1);
              loc_on_this_node.push_back(pos2 - y);
            }
            else
            {
              int cut = end - z + 1;
              new_seq = "N";
              read_on_this_node.push_back(-1);
              loc_on_this_node.push_back(-cut);
            }

            // add node
            int new_node_id = p_seq.size();
            p_seq.push_back(new_seq);
            p_read_on_node.push_back(read_on_this_node);
            p_pos_on_node.push_back(loc_on_this_node);

            BoostSTRVertex new_node;
            if(vertex.count(new_node_id) == 0)
            {
              STRVertexType node(new_node_id);
              new_node = boost::add_vertex(node, *p_graph_);
              vertex[new_node_id] = new_node;
            }
            else
            {
              new_node = vertex[new_node_id];
            }

            BoostSTRVertex v_source = vertex[node1];
            std::pair<BoostSTREdge, bool> e_search1 = boost::add_edge(v_source, new_node, *p_graph_);
            if(!e_search1.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";

            BoostSTRVertex v_target = vertex[node2];
            std::pair<BoostSTREdge, bool> e_search2 = boost::add_edge(new_node, v_target, *p_graph_);
            if(!e_search2.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
          }
        }
      }
      // after last
      pos1 = it1->first;
      node1 = it1->second;
      int e = pos1 + f_match[ref_id][node1]-1;
      while((e < coverage_2D[ref_id].size()) && coverage_2D[ref_id][e]) e++;
      if(e != (pos1 + f_match[ref_id][node1]-1))
      {
        int last_read_id = p_read_on_node[node1][0];
        int s = pos1 + f_match[ref_id][node1] - p_read_length[last_read_id] - 1;
        int len = e - s;
        std::string new_seq = seq.substr(s, len);
        IntegerVector read_on_this_node, loc_on_this_node;
        read_on_this_node.push_back(last_read_id);
        loc_on_this_node.push_back(1);

        int new_node_id = p_seq.size();
        p_seq.push_back(new_seq);
        p_read_on_node.push_back(read_on_this_node);
        p_pos_on_node.push_back(loc_on_this_node);

        BoostSTRVertex new_node;
        if(vertex.count(new_node_id) == 0)
        {
          STRVertexType node(new_node_id);
          new_node = boost::add_vertex(node, *p_graph_);
          vertex[new_node_id] = new_node;
        }
        else
        {
          new_node = vertex[new_node_id];
        }

        BoostSTRVertex v_source = vertex[node1];
        std::pair<BoostSTREdge, bool> e_search1 = boost::add_edge(v_source, new_node, *p_graph_);
        if(!e_search1.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";

      }
    }
  }
  /*********/
  t = clock() - t;
  std::cout << "merge graph takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /*********/

  delete ds;
}

bool StrGraph::loadAnchor(
  std::string& tbloutname,
  int MAX_EXTEND_LENGTH,
  anchors& anchor_set)
{
  std::ifstream fin;
  std::string cmname = "";
  fin.open(tbloutname.c_str());
  if (!fin.is_open())	// fail to open
  {
    std::cout << "Error: unable to open " << tbloutname << "\n";
    return false;
  }
  else				// succeed to open
  {
    std::string line;
    while(std::getline(fin, line))
    {
      if(line.find("Query file") != std::string::npos)
      {
        int p = line.find(":");
        while(line[++p] == ' ');
        cmname = line.substr(p);
      }
    }
  }

  fin.clear();
  fin.close();
  Integer2Integer clen_of_rfam;
  fin.open(cmname.c_str());
  if (!fin.is_open())	// fail to open
  {
    std::cout << "Error: unable to open " << cmname << "\n";
    return false;
  }
  else
  {
    std::string line;
    int clen = 0, rf = 0;
    while(std::getline(fin, line))
    {
      if(line.find("ACC") != std::string::npos)
      {
        std::string temp1 = line.substr(line.find("RF")+2);
        rf = stoi(temp1);
      }
      if(line.find("CLEN") != std::string::npos)
      {
        std::string temp2 = line.substr(line.find(" ")+1);
        clen = std::stoi(temp2);
        if(clen_of_rfam.count(rf) == 0) clen_of_rfam[rf] = clen;
      }
    }
  }
  fin.clear();
  fin.close();


  std::unordered_map<int, anchors> anchor_per_family;
  fin.open(tbloutname.c_str());
  if (!fin.is_open())	// fail to open
  {
    std::cout << "Error: unable to open " << tbloutname << "\n";
    return false;
  }
  else				// succeed to open
  {
    std::string line;
    while(std::getline(fin, line))
    {
      if(line[0] != '#')
      {
        struct anchor temp;
        if(parseAnchorLine(line, temp, clen_of_rfam, MAX_EXTEND_LENGTH))
        {

          anchor_set.push_back(temp);
          // if(ff_aggressive)
          // {
          //   int RF = temp.RF;
          //   if(anchor_per_family.count(RF) > 0) anchor_per_family[RF].push_back(temp);
          //   else                                anchor_per_family[RF] = anchors{temp};
          // }
          // else
          // {
          //   anchor_set.push_back(temp);
          // }
        }
      }
    }
  }
  fin.clear();
  fin.close();

  // if(ff_aggressive) sortAnchor(anchor_per_family, anchor_set);
  return true;
}

bool StrGraph::parseAnchorLine(
  const std::string& line,
  struct anchor& new_anchor,
  Integer2Integer& clen_of_rfam,
  int MAX_EXTEND_LENGTH)
{
  boost::char_separator<char> sep(" ");
  boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
  boost::tokenizer<boost::char_separator<char> >::iterator it;
  /* CM search */
  it = tokens.begin();  // query name   : 595,11855,1,2003,6220,1948,246:593,594
  int node = std::stoi((*it));
  /**********************/
  BoostSTRVertex this_node = vertex[node];
  (*p_graph_)[this_node].ff_seed = true;         // p_crossed[node] = true;
  (*p_graph_)[this_node].ff_delete = false;      // p_crossed[node] = true;
  /**********************/
  int length = p_seq[node].size();
  ++it;                 // accession    : -
  ++it;                 // target name  : 5_8S_rRNA
  ++it;                 // accession    : RF00002
  int RF = std::stoi((*it).substr(2));
  if(clen_of_rfam.count(RF) > 0)
  {
    int CLEN = clen_of_rfam[RF];
    int margin = CLEN / 10 + 1;
    ++it;                 // mdl          : cm
    ++it;                 // mdl from     : 38
    int mdl_start = std::stoi((*it));
    ++it;                 // mdl to       : 154
    int mdl_end = std::stoi((*it));
    ++it;                 // seq from     : 2
    int seq_start = std::stoi((*it));
    ++it;                 // seq to       : 113
    int seq_end = std::stoi((*it));
    ++it;                 // strand       : +
    ++it;                 // trunc        : no
    ++it;                 // pass         : 2
    ++it;                 // gc           : 0.48
    ++it;                 // bias         : 0.0
    ++it;                 // score        : 7.1
    ++it;                 // E-value      : 0.13
    double e_value = parseEvalue((*it));

    int left, right;
    if(seq_start < seq_end) // + strand
    {
      left = mdl_start + margin - seq_start;
      right = seq_end + CLEN + margin - mdl_end - length;
    }
    else  // - strand
    {
      left = CLEN - mdl_end + margin - seq_end;
      right = mdl_start + margin + seq_start - length;
    }

    if(MAX_EXTEND_LENGTH >= 0)
    {
      if(left   > MAX_EXTEND_LENGTH)  left = MAX_EXTEND_LENGTH;
      if(right  > MAX_EXTEND_LENGTH) right = MAX_EXTEND_LENGTH;
    }

    new_anchor.path = IntegerList{node};
    new_anchor.M_upstream = left;
    new_anchor.M_downstream = right;
    new_anchor.RF = RF;
    new_anchor.E = e_value;
    return true;
  }
  return false;
}

bool StrGraph::DirectedDFS(
  std::string& tbloutname,
  int MAX_EXTEND_LENGTH,
  std::string& path_name,
  bool aggressive,
  int batchsize)
{
  /* load Anchor */
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  anchors anchor_set;
  if(!loadAnchor(tbloutname, MAX_EXTEND_LENGTH, anchor_set)) return false;
  /** timer **/
  t = clock() - t;
  std::cout << "load anchor takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  t = clock();
  /***********/

  /* clear output */
  std::ofstream path_out;
  path_out.open(path_name.c_str());
  path_out.clear();
  path_out.close();
  p_node_id = 0;

  // temp output
  std::string temp_path_name;
  temp_path_name = "temp." + path_name;
  path_out.open(temp_path_name.c_str());
  path_out.clear();
  path_out.close();

  for(auto it = anchor_set.begin(); it != anchor_set.end(); ++it)
  {
    DirectedDFS((*it), path_name, temp_path_name, aggressive, batchsize);

    int anchor_id = (*it).path.front();
    BoostSTRVertex anchor_node = vertex[anchor_id];
    if(p_mask_edge_for_anchor.count( anchor_id ) > 0) p_mask_edge_for_anchor.erase(anchor_id);
    (*p_graph_)[anchor_node].ff_delete = true;        // this anchor has been extened
  }

  std::string cmd_cdhit = DIR_CDHIT;
  cmd_cdhit += " -M 0 -c 0.99 -T " + CPU + " -o out." + temp_path_name + " -i " + temp_path_name;
  system(cmd_cdhit.c_str());
  std::string cmd_cat = "cat out." + temp_path_name + " >> " + path_name;
  system(cmd_cat.c_str());
  std::string cmd_rm = "rm -f out.* " + temp_path_name;
  system(cmd_rm.c_str());

  /** timer **/
  t = clock() - t;
  std::cout << "extend anchor takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  return true;
}

void StrGraph::DirectedDFS(struct anchor& anchor, std::string& path_name, std::string& temp_path_name, bool aggressive, int batchsize)
{
  int anchor_id = anchor.path.front();
  std::stack<struct anchor> cheking_stacks, half_done_stacks;
  cheking_stacks.push(anchor);
  // go up first
  while(!cheking_stacks.empty())
  {
    struct anchor top_anchor = cheking_stacks.top();
    cheking_stacks.pop();
    if(top_anchor.M_upstream > 0) // need go deeper
    {
      IntegerList path = top_anchor.path;
      int first_node_id = path.front();
      BoostSTRVertex first_node = vertex[first_node_id];
      if(boost::in_degree(first_node, *p_graph_) == 0)
      {
        half_done_stacks.push(top_anchor);
      }
      else
      {
        bool ff_maskall = true;
        auto it_inegde = boost::in_edges(first_node, *p_graph_).first;
        for(; it_inegde != boost::in_edges(first_node, *p_graph_).second; ++it_inegde)
        {
          int edge_id = (*p_graph_)[*it_inegde].sid_;
          if( (p_mask_edge_for_anchor.count(anchor_id) == 0) || (p_mask_edge_for_anchor.count(anchor_id) > 0 && p_mask_edge_for_anchor[anchor_id].count(edge_id) == 0 ))
          {
            ff_maskall = false;
            BoostSTRVertex source_vertex = boost::source(*it_inegde, *p_graph_);
            int new_head = (*p_graph_)[source_vertex].rid_;
            if( !formCycle(path, new_head) )
            {
              /* mask edge */
              if( (aggressive) && ((*p_graph_)[source_vertex].ff_seed) && (!(*p_graph_)[source_vertex].ff_delete) )
              {
                if(p_mask_edge_for_anchor.count(new_head) == 0) p_mask_edge_for_anchor[new_head] = IntegerSet({edge_id});
                else                                            p_mask_edge_for_anchor[new_head].insert(edge_id);
              }
              /*************/
              path.push_front(new_head);
              int last_read_length;
              IntegerVector reads = p_read_on_node[new_head];
              int last_read_id = *(reads.end()-1);
              last_read_length = p_read_length[last_read_id];
              // update anchor stack
              struct anchor new_anchor;
              new_anchor.path = path;
              new_anchor.M_upstream = top_anchor.M_upstream - p_seq[new_head].size() + last_read_length;
              new_anchor.M_downstream = top_anchor.M_downstream;
              new_anchor.RF = top_anchor.RF;
              new_anchor.E = top_anchor.E;
              cheking_stacks.push(new_anchor);
              path.pop_front();
            }
          }
        }

        if(ff_maskall)  half_done_stacks.push(top_anchor);

      }
    }
    else                          // no need to go deeper
    {
      half_done_stacks.push(top_anchor);
    }
  }
  // go down
  while(!half_done_stacks.empty())
  {
    struct anchor top_anchor = half_done_stacks.top();
    half_done_stacks.pop();
    if(top_anchor.M_downstream > 0) // can go deeper
    {
      IntegerList path = top_anchor.path;
      int last_node_id = path.back();
      BoostSTRVertex last_node = vertex[last_node_id];
      if(boost::out_degree(last_node, *p_graph_) == 0)
      {
        savePath(path_name, temp_path_name, top_anchor, batchsize);
      }
      else
      {
        bool ff_maskall = true;
        auto it_outegde = boost::out_edges(last_node, *p_graph_).first;
        for(; it_outegde != boost::out_edges(last_node, *p_graph_).second; ++it_outegde)
        {
          int edge_id = (*p_graph_)[*it_outegde].sid_;
          if( (p_mask_edge_for_anchor.count(anchor_id) == 0) || (p_mask_edge_for_anchor.count(anchor_id) > 0 && p_mask_edge_for_anchor[anchor_id].count(edge_id) == 0 )) // this edge is not masked
          {
            ff_maskall = false;
            BoostSTRVertex target_vertex = boost::target(*it_outegde, *p_graph_);
            int new_tail = (*p_graph_)[target_vertex].rid_;
            if( !formCycle(path, new_tail) )
            {
              /* mask edge */
              if( (aggressive) && ((*p_graph_)[target_vertex].ff_seed) && (!(*p_graph_)[target_vertex].ff_delete) )
              {
                if(p_mask_edge_for_anchor.count(new_tail) == 0) p_mask_edge_for_anchor[new_tail] = IntegerSet({edge_id});
                else                                            p_mask_edge_for_anchor[new_tail].insert(edge_id);
              }
              /*************/
              path.push_back(new_tail);
              int first_read_length;
              IntegerVector reads = p_read_on_node[new_tail];
              int first_read_id = reads[0];
              first_read_length = p_read_length[first_read_id];

              // update anchor stack
              struct anchor new_anchor;
              new_anchor.path = path;
              new_anchor.M_upstream = top_anchor.M_upstream;
              new_anchor.M_downstream = top_anchor.M_downstream - p_seq[new_tail].size() + first_read_length;
              new_anchor.RF = top_anchor.RF;
              new_anchor.E = top_anchor.E;
              half_done_stacks.push(new_anchor);
              path.pop_back();
            }
          }
        }
        if(ff_maskall)  savePath(path_name, temp_path_name, top_anchor, batchsize);
      }
    }
    else
    {
      savePath(path_name, temp_path_name, top_anchor, batchsize);
    }
  }
}

void StrGraph::savePath(std::string& path_name, std::string& temp_path_name, const struct anchor& anchor, int batchsize)
{
  ++p_node_id;
  std::ofstream path_out;
  path_out.open(temp_path_name.c_str(),std::ios::in|std::ios::out|std::ios::app);

  IntegerList path = anchor.path;
  std::string seq;
  path_out << ">" << "RF" << anchor.RF << "_" << p_node_id;
  for(auto it = path.begin(); it != path.end(); ++it) path_out << "," << (*it);
  path_out << " ";
  // eg: path_i,node_x1,node_x2,...,node_xn r_y1,begin,length,r_y2,begin,length,..,r_ym,begin,length,
  // all integer, don't cut here
  int M = 0;

  for(auto it = path.begin(); it != path.end(); ++it)
  {
    int edge = *it;
    // /*********************/
    // std::cout << edge << "\n";
    // /*********************/
    // assembly sequence
    std::string full = p_seq[edge];
    int first_read = p_read_on_node[edge][0];
    if(it != path.begin())
    {
      if(first_read < 0)  // merged node
      {
        seq = seq.substr(0, seq.size() + p_pos_on_node[edge][0]);
      }
      else
      {
        seq += full.substr(p_read_length[first_read]);
      }
    }
    else                   seq += full;
    // print reads on path
    IntegerVector reads = p_read_on_node[edge];
    IntegerVector locus = p_pos_on_node[edge];
    if(it != path.begin())
    {
      for(int i = 1; i < reads.size(); i++)
      {
        path_out  << reads[i]                     << ","
                  << std::to_string(locus[i] + M) << ",";
        if(reads[i] < 0) path_out << "-1,";
        else             path_out << p_read_length[reads[i]]  << ",";
      }
    }
    else
    {
      for(int i = 0; i < reads.size(); i++)
      {
        path_out  << reads[i]                     << ","
                  << std::to_string(locus[i] + M) << ",";
        if(reads[i] < 0) path_out << "-1,";
        else             path_out << p_read_length[reads[i]]  << ",";
      }
    }
    if(full[0] != 'N') M += full.size() - p_read_length[*(reads.end()-1)];
  }
  path_out << "\n" << seq << "\n";
  path_out.close();

  if(p_node_id > 0 && (p_node_id % batchsize == 0))
  {
    std::string cmd_cdhit = DIR_CDHIT;
    cmd_cdhit += " -M 0 -c 0.99 -T " + CPU + " -o out." + temp_path_name + " -i " + temp_path_name;
    system(cmd_cdhit.c_str());
    std::string cmd_cat = "cat out." + temp_path_name + " >> " + path_name;
    system(cmd_cat.c_str());
    std::string cmd_rm = "rm -f out.* " + temp_path_name;
    system(cmd_rm.c_str());
  }
}

void StrGraph::writeMergedGraph(std::string& filename)
{
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  std::cout << "start write..\n";
  std::string outname = filename.substr(0, filename.find_last_of('.')) + ".merged.fq";
  std::ofstream fout(outname.c_str());
  for(int node_i = 0; node_i < p_seq.size(); node_i++)
  {
    fout << ">" << node_i << ",";
    for(int r = 0; r < p_read_on_node[node_i].size(); r++)
    {
      fout << p_read_on_node[node_i][r] << ',' << p_pos_on_node[node_i][r] << ',';
      if(p_read_on_node[node_i][r] < 0)        fout << "-1";
      else                                     fout << p_read_length[p_read_on_node[node_i][r]];
      if(r+1 < p_read_on_node[node_i].size())  fout << ',';
    }

    BoostSTRVertex v = vertex[node_i];
    auto it_v = boost::adjacent_vertices(v, *p_graph_).first;
    if(it_v != boost::adjacent_vertices(v, *p_graph_).second)
    {
      fout << ':';
      for(; it_v != boost::adjacent_vertices(v, *p_graph_).second; it_v++)
      {
        fout << (*p_graph_)[*it_v].rid_;
        it_v++;
        if(it_v != boost::adjacent_vertices(v, *p_graph_).second) fout << ',';
        it_v--;
      }
    }
    fout << "\n" << p_seq[node_i] << "\n";
  }
  /** timer **/
  t = clock() - t;
  std::cout << "write graph takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  fout.close();
}

/***************************** for overlap graph ******************************/
bool StrGraph::readAsqgFile(std::string& filename)
{
  std::ifstream fin;
  // fin.open(config.IDIR_GRAPH.c_str());	// open file
  fin.open(filename.c_str());	// open file
  if (!fin.is_open())	    // fail to open
  {
    std::cout << "Error: " << filename << " doesnot exist!\n";
    return false;
  }
  else				// succeed to open
  {
    std::string line;
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    /** timer **/
    clock_t t;
    t = clock();
    /***********/
    while(std::getline(fin, line))
    {
      if(line[0] == 'V')
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();            // cont[0], VT
        ++it;                           // cont[1], read name
        int id = std::stoi(*it);
        ++it;                           // cont[2], sequence
        p_seq.push_back(*it);
      }
      if(line[0] == 'E')      break;
    }
    /** timer **/
    t = clock() - t;
    std::cout << "VT takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    t = clock();
    /***********/

    p_order = p_seq.size();
    p_traversed = BooleanVector(2*p_order, false);
    p_crossed   = BooleanVector(2*p_order, false);
    /** timer **/
    t = clock() - t;
    std::cout << "resize takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    t = clock();
    /***********/
    do
    {
      if(line[0] == 'E')  // edge, sperate by '\t' and ' '
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();            // cont[0], ED
        ++it;                           // cont[1], read1 name
        int source_id = std::stoi(*it);
        ++it;                           // cont[2], read2 name
        int target_id = std::stoi(*it);
        ++it;                           // cont[3], A
        int a = std::stoi(*it);
        ++it;                           // cont[4], B
        int b = std::stoi(*it);
        ++it;                           // cont[5], l1
        int l1 = std::stoi(*it);
        ++it;                           // cont[6], C
        int c = std::stoi(*it);
        ++it;                           // cont[7], D
        int d = std::stoi(*it);
        ++it;                           // cont[8], l2
        int l2 = std::stoi(*it);
        ++it;                           // cont[9], "1" == RC
        bool ff_RC = *it == "1";

        /* add or locate the two vertices in the graph */
        BoostSTRVertex v_source, v_RC_source, v_target, v_RC_target;
        if(vertex.count(source_id) == 0)
        {
          STRVertexType node(source_id);
          node.len_ = l1;
          v_source = boost::add_vertex(node, *p_graph_);
          vertex[source_id] = v_source;
        }
        else
        {
          v_source = vertex[source_id];
        }

        if(vertex.count(source_id + p_order) == 0)
        {
          STRVertexType node(source_id + p_order);
          node.len_ = l1;
          v_RC_source = boost::add_vertex(node, *p_graph_);
          vertex[source_id + p_order] = v_RC_source;
        }
        else
        {
          v_RC_source = vertex[source_id + p_order];
        }

        if(vertex.count(target_id) == 0)
        {
          STRVertexType node(target_id);
          node.len_ = l2;
          v_target = boost::add_vertex(node, *p_graph_);
          vertex[target_id] = v_target;
        }
        else
        {
          v_target = vertex[target_id];
        }

        if(vertex.count(target_id + p_order) == 0)
        {
          STRVertexType node(target_id + p_order);
          node.len_ = l2;
          v_RC_target = boost::add_vertex(node, *p_graph_);
          vertex[target_id + p_order] = v_RC_target;
        }
        else
        {
          v_RC_target = vertex[target_id + p_order];
        }

        if(source_id != target_id)
        {
          /* add edge */
          if(ff_RC)         // r2 is reversed
          {
            if(a != 0)      // r1->r2
            {
              /* r1 -> r2' */
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_RC_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = a;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              /* r2 -> r1' */
              e_search = boost::add_edge(v_target, v_RC_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = c;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
            }
            else            // r2->r1
            {
              /* r2' -> r1 */
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_RC_target, v_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l2-1-d;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              /* r1' -> r2 */
              e_search = boost::add_edge(v_RC_source, v_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l1-1-b;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
            }
          }
          else
          {
            if(a != 0)     // r1->r2
            {
              /* r1 -> r2 */
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = a;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              /* r2' -> r1' */
              e_search = boost::add_edge(v_RC_target, v_RC_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l2-1-d;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
            }
            else            // r2->r1
            {
              /* r2 -> r1 */
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_target, v_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = c;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              /* r1' -> r2' */
              e_search = boost::add_edge(v_RC_source, v_RC_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l1-1-b;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
            }
          }
        }
      }
    }
    while(std::getline(fin, line));
    /** timer **/
    t = clock() - t;
    std::cout << "ED takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    /***********/

  }
  fin.clear();
  fin.close();
  return true;
}

void StrGraph::CondenseGraph()
{
  /** timer **/
  clock_t t;
  /***********/
  /* find seed vertex */
  t = clock();
  int total_node = 0;
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++it)
  {
    ++total_node;
    if(boost::in_degree(*it, *p_graph_) <= 0)
    { // orphant read is not seed
      if(boost::out_degree(*it, *p_graph_) > 0)   (*p_graph_)[*it].ff_seed = true;
    }
    else
    {
      if(boost::out_degree(*it, *p_graph_) != 1)  (*p_graph_)[*it].ff_seed = true;
      if(boost::in_degree(*it, *p_graph_) > 1)    (*p_graph_)[*it].ff_seed = true;
    }
  }
  t = clock() - t;
  std::cout << total_node << " vertices in total\n";
  std::cout << "find seed takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";


  /* for each component */
  std::list<cycle_s> to_add_cycle;
  int cnt_cycle_node = 0;
  double condeseSE = 0, condenseCycle = 0, findComponent = 0;
  BooleanVector connected = BooleanVector(2*p_order, false);
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++it)
  {
    if(! p_crossed[ (*p_graph_)[*it].rid_ ]) // untraversed vertex, new component
    {
      /* find connected component and source_edges */
      t = clock();
      std::stack<BoostSTRVertex> to_visit;
      to_visit.push(*it);
      std::list<BoostSTRVertex> this_component;
      int cnt_vertex = 0;
      while(!to_visit.empty())
      {
        BoostSTRVertex top_vertex = to_visit.top();
        to_visit.pop();
        int rid = (*p_graph_)[top_vertex].rid_;
        if(!connected[rid])
        {
          connected[rid] = true;
          this_component.push_back(top_vertex);
          ++cnt_vertex;
          p_crossed[rid] = true;
          // if(rid < p_order) p_crossed[rid+p_order] = true;
          // else              p_crossed[rid-p_order] = true;
          for(auto it_v = boost::adjacent_vertices(top_vertex, *p_graph_).first; it_v != boost::adjacent_vertices(top_vertex, *p_graph_).second; ++it_v)
          {
            if(!connected[ (*p_graph_)[*it_v].rid_ ]) to_visit.push(*it_v);
          }
          for(auto it_v = boost::inv_adjacent_vertices(top_vertex, *p_graph_).first; it_v != boost::inv_adjacent_vertices(top_vertex, *p_graph_).second; ++it_v)
          {
            if(!connected[ (*p_graph_)[*it_v].rid_ ]) to_visit.push(*it_v);
          }
        }
      }
      t = clock() - t;
      findComponent += ((double)t)/CLOCKS_PER_SEC;

      /* condense this component */
      t = clock();
      bool ff_singleCycle = true;
      for(auto it_v = this_component.begin(); it_v != this_component.end(); ++it_v)
      {
        if((*p_graph_)[*it_v].IsSeed() && (!p_traversed[ (*p_graph_)[*it_v].rid_ ]))
        {
          ff_singleCycle = false;
          condense(*it_v, to_add_cycle);
        }
      }
      t = clock() - t;
      condeseSE += ((double)t)/CLOCKS_PER_SEC;

      t = clock();
      if(ff_singleCycle)  // if this component is single loop
      {
        cnt_cycle_node += cnt_vertex;
        auto it_v = this_component.begin();
        (*p_graph_)[*it_v].ff_seed = true;
        condense(*it_v, to_add_cycle);
      }
      t = clock() - t;
      condenseCycle += ((double)t)/CLOCKS_PER_SEC;
    }
  }
  std::cout << "find components takes " << findComponent << " s\n";
  std::cout << "consense SE takes " << condeseSE << " s\n";
  std::cout << "consense cycle takes " << condenseCycle << " s\n";
  std::cout << cnt_cycle_node << " cycle nodes\n";
  /* modify graph */
  t = clock();
  std::list<BoostSTREdge> to_delete_edge;
  for(auto it_e = boost::edges(*p_graph_).first; it_e != boost::edges(*p_graph_).second; ++it_e)
  {
    if( (*p_graph_)[*it_e].ff_delete ) to_delete_edge.push_back(*it_e);
  }
  std::list<BoostSTRVertex> to_delete_vertex;
  for(auto it_v = boost::vertices(*p_graph_).first; it_v != boost::vertices(*p_graph_).second; ++it_v)
  {
    if( (*p_graph_)[*it_v].ff_delete ) to_delete_vertex.push_back(*it_v);
  }
  for(auto it = to_delete_edge.begin(); it != to_delete_edge.end(); ++it)     boost::remove_edge(*it, *p_graph_);
  for(auto it = to_delete_vertex.begin(); it != to_delete_vertex.end(); ++it) boost::remove_vertex(*it, *p_graph_);

  p_node_id = 0;
  for(auto it = to_add_cycle.begin(); it != to_add_cycle.end(); ++it)
  {
    BoostSTRVertex s = it->head;
    BoostSTRVertex t = it->tail;
    std::pair<BoostSTREdge, bool> edge_new = add_edge(s, t, *p_graph_);
    if(edge_new.second)
    {
      (*p_graph_)[edge_new.first].SetCondensedTag(true);
      (*p_graph_)[edge_new.first].path_info_ = it->path;
      (*p_graph_)[edge_new.first].ff_cycle   = it->ff_cycle;
      (*p_graph_)[edge_new.first].sid_       = p_node_id;
      ++p_node_id;
    }
  }

  t = clock() - t;
  std::cout << "modify graph takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
}

void StrGraph::condense(const BoostSTRVertex seed,
  std::list<cycle_s>& to_add_cycle)
{
  for(auto it_e = boost::out_edges(seed, *p_graph_).first; it_e != boost::out_edges(seed, *p_graph_).second; ++it_e)
    condense(*it_e, to_add_cycle);
}

void StrGraph::condense(const BoostSTREdge source_edge,
  std::list<cycle_s>& to_add_cycle)
{
  BoostSTREdge current_edge = source_edge;
  BoostSTRVertex head = boost::source(current_edge, *p_graph_);
  BoostSTRVertex tail = head;
  IntegerVector path_info;
  path_info.push_back((*p_graph_)[head].rid_);
  path_info.push_back((*p_graph_)[head].len_);
  while(1)
  {
    p_traversed[ (*p_graph_)[tail].rid_ ] = true;
    (*p_graph_)[current_edge].ff_delete = true;             //boost::remove_edge(current_edge, *p_graph_);
    if(tail != head)  (*p_graph_)[tail].ff_delete = true;   // boost::remove_vertex(to_delete, *p_graph_);

    tail = boost::target(current_edge, *p_graph_);          // define the new tail vertex
    path_info.push_back((*p_graph_)[current_edge].p_A);
    path_info.push_back((*p_graph_)[tail].rid_);
    path_info.push_back((*p_graph_)[tail].len_);

    if((*p_graph_)[tail].IsSeed())  break;
    else                            current_edge = *(boost::out_edges(tail, *p_graph_)).first;
  }
  path_info.push_back(-1);
  if( boost::out_degree(tail, *p_graph_) == 0)  p_traversed[ (*p_graph_)[tail].rid_ ] = true;

  cycle_s new_cycle;
  new_cycle.head = head;
  new_cycle.tail = tail;
  new_cycle.path = path_info;
  new_cycle.ff_cycle = ((*p_graph_)[head].rid_ == (*p_graph_)[tail].rid_ );
  to_add_cycle.push_back(new_cycle);
}

void StrGraph::writeGraph(std::string& filename)
{
  /** timer **/
  clock_t t;
  t = clock();
  int cnt = 0;
  /***********/
  std::string newname = filename;
  std::string outname = newname.substr(0,newname.find_last_of('.')) + ".StringGraph.fq";
  std::ofstream fout(outname.c_str());
  for(auto it_e = boost::edges(*p_graph_).first; it_e != boost::edges(*p_graph_).second; ++it_e)
  {
    if((*p_graph_)[*it_e].IsCondensed())
    {
      ++cnt;
      std::string sequence = "";
      IntegerVector path_info = (*p_graph_)[*it_e].path_info_;
      fout << ">" << (*p_graph_)[*it_e].sid_ << ",";
      int p_before = 1;
      for(int i = 0; i < path_info.size(); i += 3)
      {
        int r_id = path_info[i];
        int r_len = path_info[i+1];
        int a = path_info[i+2];
        std::string temp = (r_id < p_order) ? p_seq[r_id] : RC(r_id-p_order);
        if(a > 0)   sequence += temp.substr(0, a);
        else        sequence += temp;

        fout << ((r_id < p_order) ? r_id : r_id-p_order) << ","
             << p_before << ","
             << r_len;
        if(i + 3 < path_info.size())  fout << ",";
        p_before = a + p_before;

        if(i + 6 >= path_info.size() && (*p_graph_)[*it_e].IsCycle())
        {
          i = -3;
          (*p_graph_)[*it_e].ff_cycle = false;
        }
      }

      BoostSTRVertex t = boost::target(*it_e, *p_graph_);
      auto it_out_e = boost::out_edges(t, *p_graph_).first;
      if(it_out_e != boost::out_edges(t, *p_graph_).second)
      {
        fout << ":" << (*p_graph_)[*it_out_e].sid_;
        ++it_out_e;
        while(it_out_e != boost::out_edges(t, *p_graph_).second)
        {
          fout << "," << (*p_graph_)[*it_out_e].sid_;
          ++it_out_e;
        }
      }
      fout << "\n";
      fout << sequence << "\n";
    }
  }
  /** timer **/
  t = clock() - t;
  std::cout << "write "<< cnt << " edges " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  t = clock();
  cnt = 0;
  /***********/

  /* take care of orphant reads */
  for(int i = 0; i < p_order; ++i)
  {
    if(!p_traversed[i] && !p_traversed[i+p_order])
    {
      ++cnt;
      fout << ">" << p_node_id << "," << i << ",1," << p_seq[i].size() << "\n";
      fout << p_seq[i] << "\n";
      ++p_node_id;
    }
  }

  /** timer **/
  t = clock() - t;
  std::cout << "write "<< cnt << " orphants " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  fout.close();
}

/*************************** auxiliary function *******************************/
StrGraph::StrGraph()
{
  p_graph_ = new BoostSTRGraph();
}

StrGraph::~StrGraph()
{
  delete p_graph_;
}

void StrGraph::showGraph()
{
  std::cout << "seq:\n";
  for(int i = 0; i < p_seq.size(); ++i)     std::cout << i << " : " << p_seq[i] << "\n";
  std::cout << "\n";

  std::cout << "vertices:\n";
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++it)
  {
    std::cout << "Vertex " << (*p_graph_)[*it].rid_;
    if(!ff_stringGraph) std::cout << ", len = " << (*p_graph_)[*it].len_;
    std::cout << "\n";
  }
  std::cout << "\n";

  std::cout << "edges:\n";
  for(auto it_e = boost::edges(*p_graph_).first; it_e != boost::edges(*p_graph_).second; ++it_e)
  {
    BoostSTRVertex s = boost::source(*it_e, *p_graph_);
    BoostSTRVertex t = boost::target(*it_e, *p_graph_);
    std::cout << "Edge: ("<< (*p_graph_)[s].rid_ << ", " << (*p_graph_)[t].rid_ << ")";
    if(!ff_stringGraph)  std::cout << ": " << (*p_graph_)[*it_e].p_A;
    std::cout << "\n";
  }
  std::cout << "\n";

  if(ff_stringGraph)
  {
    std::cout << "read_length:\n";
    for(auto x : p_read_length) std::cout << x.first << ": " << x.second << "\n";
    std::cout << "\n";

    std::cout << "edge_info:\n";
    for(int i = 0; i < p_order; ++i)
    {
      for(auto x : p_read_on_node[i]) std::cout << x << " ";
      std::cout << "\n";
      for(auto x : p_pos_on_node[i])  std::cout << x << " ";
      std::cout << "\n";
    }
    std::cout << "\n";
  }
}

double parseEvalue(const std::string& e)
{
  return (e == "0") ? 999999 : -log10(std::stod(e));
}


bool formCycle(IntegerList & path, int new_tail)
{
  for(auto it = path.begin(); it != path.end(); ++it)
  {
    if((*it) == new_tail)  return true;
  }
  return false;
}

std::string StrGraph::RC(int i)
{
  std::string res="";
  for(int j=0;j<p_seq[i].size();++j)
  {
    if(p_seq[i][j]=='A' || p_seq[i][j]=='a')
      res = "T" + res;
    else if (p_seq[i][j]=='T' || p_seq[i][j]=='t')
      res = "A" + res;
    else if (p_seq[i][j]=='G' || p_seq[i][j]=='g')
      res = "C" + res;
    else if (p_seq[i][j]=='C' || p_seq[i][j]=='c')
      res = "G" + res;
  }
  return res;
}

bool compare_nocase (const struct anchor& first, const struct anchor& second)
// true if the first item goes before the second
{
  // if(first.path.front() == second.path.front()) return (first.E < second.E);
  return (first.E > second.E);
}


void sortAnchor(std::unordered_map<int, anchors>& anchorperf, anchors& anchor_set)
{
  for(auto x: anchorperf)
  {
    x.second.sort(compare_nocase);
    anchor_set.insert(anchor_set.end(), x.second.begin(), x.second.end());
    struct anchor tt;
    tt.path = IntegerList{-1};
    anchor_set.push_back(tt);
  }
}
