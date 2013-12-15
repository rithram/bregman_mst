#ifndef MINIMUM_SPANNING_TREE_HPP_
#define MINIMUM_SPANNING_TREE_HPP_

#include "data.hpp"
#include "union_find.hpp"

namespace bmst {

  // An edge of the MST
  class Edge {
  
  public:
    
    size_t u;
    size_t v;
    double weight;
    
    Edge(size_t u_in, size_t v_in, double weight_in)
      :
    u(u_in), v(v_in), weight(weight_in)
    {}
      
    Edge()
      :
    u(-1), v(-1), weight(-DBL_MAX)
    {}
  
  }; // class Edge
  
  struct EdgeSorterStruct
  {
    bool operator()(const Edge& A, const Edge& B)
    {
      return (A.weight < B.weight);
    }
  } EdgeSorter;


  template<typename T, class EdgePolicy, class TTreeType>
  class MinimumSpanningTree {

  private:

    Table<T> data_;

    TTreeType* tree_;
    
    std::vector<size_t> old_from_new_;
    
    std::vector<Edge> edge_list_;
    
    UnionFind components_;
    
    std::vector<Edge> nearest_neighbors_;
    std::vector<double> candidate_dists_;
    
    // functions //
    
    void SearchTree_(TTreeType* query_node, TTreeType* reference_node);
    
    void NaiveBoruvka_(std::vector<std::vector<double> >& edge_weights);

    void AddEdges_();
    
    void SearchTree_(const Point<T>& q, size_t q_index, size_t root_q, TTreeType* node);
    
    void ResetTree_(TTreeType* node);

    void UpdateTree_(TTreeType* node);
    
    void ResetAll_();
    
  public:
    
    MinimumSpanningTree(Table<T>& data, int leaf_size = 1);
  
    ~MinimumSpanningTree();
  
    // Dual-tree Boruvka
    void ComputeDTB();
    
    // Naive
    // flag indicates whether we compute all edge weights and store or compute 
    // as needed
    void ComputeNaive(bool use_n2_memory);
    
    // Single-tree Boruvka, loop over all queries
    void ComputeSTB();
  
    // This will put things back in terms of the original indexing  
    std::vector<Edge>& EdgeList();
  
  }; // class MST


}


#include "minimum_spanning_tree_impl.hpp"





#endif

