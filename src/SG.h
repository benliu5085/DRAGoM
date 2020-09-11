#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_
#include <boost/tokenizer.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <list>
#include <stack>
#include <time.h>
#include <math.h>       // log10

typedef std::list<int>                        IntegerList;
struct anchor
{
  IntegerList path;
  int M_upstream = -1;
  int M_downstream = -1;
  int RF = -1;
  double E = -1;
};
bool compare_nocase (const struct anchor& first, const struct anchor& second);

typedef std::list<struct anchor>              anchors;
typedef std::vector<int>                      IntegerVector;
typedef std::vector<IntegerVector>            IntegerVector2D;
typedef std::vector<bool>                     BooleanVector;
typedef std::vector<BooleanVector>            BooleanVector2D;
typedef std::vector<std::string>              StringVector;
typedef std::unordered_map<int, std::string>  Integer2String;
typedef std::unordered_map<std::string, int>  String2Integer;
typedef std::unordered_map<int, int>          Integer2Integer;
typedef std::unordered_set<int>               IntegerSet;
typedef std::unordered_map<int, IntegerSet>   Integer2IntegerSet;

/*
adjacency_list<         // a 2D structure, the first is a vertex, vertex contains its edge list
   OutEdgeList,         // container is used to represent the edge lists, vecS faster to iterate, listS faster to add
   VertexList,          // container is used to represent the outer two-dimensional container
   Directed,            // directed / undirected / (bidirectional) directed with access to both the in-edges and out-edges
   VertexProperties,    //
   EdgeProperties,      //
   GraphProperties,     //
   EdgeList>            //
*/
class STRVertexType
{
public:
  explicit STRVertexType(const int i)       { rid_ = i;       }
  inline bool IsSeed(void)                  { return ff_seed; }
  /* use for overlap graph */
  int rid_ = -1;                // the read ID for the current vertex
  int len_ = 0;                 // the length of the read
  bool ff_delete = false;       // tag to delete after condense
  bool ff_seed = false;         // tag to stop travering
};

class STREdgeType
{
public:
  inline void SetCondensedTag(const bool i) { ff_condensed = i; }
  inline bool IsCondensed(void)             { return ff_condensed; }
  inline bool IsCycle(void)                 { return ff_cycle; }
  /* use for overlap graph */
  int sid_ = -1;
  bool ff_condensed = false;             // tag to check whether the edge has been condensed
  bool ff_delete = false;                // tag to delete after condense cycle
  bool ff_cycle = false;
  IntegerVector path_info_;
  /* a vector contanins the path information of the edge with the following format:
     read_ID:read_len:p_A: read_ID:read_len:p_A ...
  */
  int p_A = -1;
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, STRVertexType, STREdgeType> BoostSTRGraph;
typedef boost::graph_traits<BoostSTRGraph>::vertex_descriptor BoostSTRVertex;
typedef boost::graph_traits<BoostSTRGraph>::edge_descriptor   BoostSTREdge;
typedef boost::graph_traits<BoostSTRGraph>::vertex_iterator   BoostSTRVertexIter;
typedef boost::graph_traits<BoostSTRGraph>::edge_iterator     BoostSTREdgeIter;

typedef struct Cycle
{
  BoostSTRVertex  head;
  BoostSTRVertex  tail;
  IntegerVector   path;
  bool            ff_cycle;
} cycle_s;


double parseEvalue(const std::string& e);
bool formCycle(IntegerList & path, int new_tail);
void sortAnchor(std::unordered_map<int, anchors>& anchorperf, anchors& anchor_set);
void concateAnchor(std::unordered_map<int, anchors>& anchorperf, anchors& anchor_set);
unsigned int goodMapping(std::string& CIAGR, int flag);
bool goodJunction(BooleanVector& coverage, int pos_i, int pos_j);
bool thisIsOrphat(std::string& header);
void accumulateDegree(const std::string& header);

class StrGraph {
private:
  bool ff_stringGraph = false;        // false if overlap graph
  bool ff_spadesGraph = false;        // true if the graph from SPAdes
  // the overlap graph
  long long int   p_node_id = 0;      // both
  int             p_order;
  BoostSTRGraph*  p_graph_;           // both
  StringVector    p_seq;              // both
  BooleanVector   p_traversed;
  BooleanVector   p_crossed;
  std::unordered_map<int, BoostSTRVertex> vertex; // both

  void markVertexAsTraversed(int v_id);
  void condense(const BoostSTRVertex seed,
    std::list<cycle_s>& to_add_cycle);
  void condense(const BoostSTREdge source_edge,
    std::list<cycle_s>& to_add_cycle);
  std::string RC(int i);

  // the string graph
  Integer2Integer     p_read_length;
  IntegerVector2D     p_read_on_node;
  IntegerVector2D     p_pos_on_node;
  Integer2IntegerSet  p_mask_edge_for_anchor;
  // BooleanVector2D masked_edge; // if enable mask
  bool loadAnchor(
    std::string& tbloutname,
    int MAX_EXTEND_LENGTH,
    anchors& anchor_set);
  bool parseAnchorLine(
      const std::string& line,
      struct anchor& new_anchor,
      Integer2Integer& clen_of_rfam,
      int MAX_EXTEND_LENGTH);
  void DirectedDFS(struct anchor& anchor, std::string& path_name, std::string& temp_path_name, bool aggressive, int batchsize);
  void savePath(std::string& path_name, std::string& temp_path_name, const struct anchor& anchor, int batchsize);

  // merge graph
public:
  StrGraph();
  ~StrGraph();
  void showGraph();

  // the overlap graph
  bool readAsqgFile(std::string& filename);
  void CondenseGraph();
  void writeGraph(std::string& filename);

  // the string graph
  bool readFqFile(std::string& filename);
  bool DirectedDFS(
    std::string& tbloutname,
    int MAX_EXTEND_LENGTH,
    std::string& path_name,
    bool aggressive,
    int batchsize);

  // merge graph
  void mergeSG(std::string& sam_name, std::string& contig_name);
  void writeMergedGraph(std::string& filename);
};
#endif
