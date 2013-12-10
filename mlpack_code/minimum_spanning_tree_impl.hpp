
#ifndef MINIMUM_SPANNING_TREE_IMPL_HPP_
#define MINIMUM_SPANNING_TREE_IMPL_HPP_

namespace bmst {

  template<typename T, class EdgePolicy, class TTreeType>
  MinimumSpanningTree<T, EdgePolicy, TTreeType>::MinimumSpanningTree(Table<T>& data)
  :
  data_(data),
  components_(data.n_points()),
  nearest_neighbors_(data.n_points()),
  candidate_dists_(data.n_points(), DBL_MAX)
  {
    
    tree_ = new TTreeType(data_, old_from_new_);
    
  }

  template<typename T, class EdgePolicy, class TTreeType>
  MinimumSpanningTree<T, EdgePolicy, TTreeType>::~MinimumSpanningTree()
  {
    delete tree_;
  }

  // Dual-tree Boruvka
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::ComputeDTB()
  {
    
    while (edge_list_.size() < data_.n_points() - 1)
    {
      
      SearchTree_(tree_, tree_);
      
      AddEdges_();
      
    }
    
  }
  
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::SearchTree_(TTreeType* query_node,
                                                                  TTreeType* reference_node)
  {
    
    if (query_node->Bound().Component() >= 0 && query_node->Bound().Component() == reference_node->Bound().Component()) 
    {
      return; // they're connected, so prune
    }
    else if (EdgePolicy::CanPrune(query_node->Bound(), reference_node->Bound()))
    {
      return; // we pruned based on bounds
    }
    else if (query_node->IsLeaf() && reference_node->IsLeaf()) {
      
      for (size_t q = query_node->Begin(); q < query_node->End(); q++)
      {
        
        const Point<T>& query = data_[q];
        size_t root_q = components_.Find(q);
        
        for (size_t r = reference_node->Begin(); reference_node->End(); r++)
        {

          const Point<T>& reference = data_[r];
          size_t root_r = components_.Find(r);
          
          double this_weight = EdgePolicy::EdgeWeight(query, reference);
          
          if (this_weight < candidate_dists_[root_q]) 
          {
            
            candidate_dists_[root_q] = this_weight;
            nearest_neighbors_[root_q] = Edge(q,r, this_weight);
            
          }
          
        } // loop over r
      } // loop over q
    } // base case
    else if (reference_node->IsLeaf())
    {
      
      SearchTree_(query_node->Left(), reference_node);
      SearchTree_(query_node->Right(), reference_node);
      
    }
    else if (query_node->IsLeaf())
    {
     
      double left_dist = EdgePolicy::EdgeWeight(query_node->Center(), reference_node->Left()->Center());
      double right_dist = EdgePolicy::EdgeWeight(query_node->Center(), reference_node->Right()->Center());
      
      if (left_dist < right_dist)
      {
        SearchTree_(query_node, reference_node->Left());
        SearchTree_(query_node, reference_node->Right());
      }
      else {
        SearchTree_(query_node, reference_node->Right());
        SearchTree_(query_node, reference_node->Left());
      }
    }
    else {
     
      double left_dist = EdgePolicy::EdgeWeight(query_node->Left()->Center(), reference_node->Left()->Center());
      double right_dist = EdgePolicy::EdgeWeight(query_node->Left()->Center(), reference_node->Right()->Center());
      
      if (left_dist < right_dist)
      {
        SearchTree_(query_node->Left(), reference_node->Left());
        SearchTree_(query_node->Left(), reference_node->Right());
      }
      else {
        SearchTree_(query_node->Left(), reference_node->Right());
        SearchTree_(query_node->Left(), reference_node->Left());
      } 

      left_dist = EdgePolicy::EdgeWeight(query_node->Right()->Center(), reference_node->Left()->Center());
      right_dist = EdgePolicy::EdgeWeight(query_node->Right()->Center(), reference_node->Right()->Center());
      
      if (left_dist < right_dist)
      {
        SearchTree_(query_node->Right(), reference_node->Left());
        SearchTree_(query_node->Right(), reference_node->Right());
      }
      else {
        SearchTree_(query_node->Right(), reference_node->Right());
        SearchTree_(query_node->Right(), reference_node->Left());
      } 

    }
    
  }
  
  // Naive
  // flag indicates whether we compute all edge weights and store or compute 
  // as needed
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::ComputeNaive(bool use_n2_memory)
  {
    
    std::vector<std::vector<double> > edge_weights;
    
    if (use_n2_memory)
    {

      edge_weights.resize(data_.n_points());  
      
      for (size_t i = 0; i < data_.n_points() - 1; i++)
      {
        
        edge_weights[i].resize(data_.n_points());
        edge_weights[i][i] = 0.0;
        
        for (size_t j = i+1; j < data_.n_points(); j++)
        {
      
          double this_weight = EdgePolicy::EdgeWeight(data_[i], data_[j]);
          edge_weights[i][j] = this_weight;
          edge_weights[j][i] = this_weight;
          
        } // for j
      } // for i
      
    }
    
    NaiveBoruvka_(edge_weights);
    
    
  }
  
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::NaiveBoruvka_(std::vector<std::vector<double> >& edge_weights)
  {
    
    bool compute_weights = (edge_weights.size() == 0);
    
    // until we have N - 1 edges
    while (edge_list_.size() < data_.n_points() - 1)
    {
      
      for (size_t i = 0; i < data_.n_points(); i++)
      {
        
        const Point<T>& point_i = data_[i];
        
        for (size_t j = 0; j < data_.n_points(); j++)
        {
        
          if (i == j) continue;
        
          const Point<T>& point_j = data_[j];
          
          if (components_.Find(i) != components_.Find(j)) 
          {
            
            size_t root_i = components_.Find(i);
            
            // they aren't in the same component
            double this_weight;
            if (compute_weights) 
            {
              this_weight = EdgePolicy::EdgeWeight(point_i, point_j);
            }
            else {
              this_weight = edge_weights[i][j];
            }
            
            if (this_weight < candidate_dists_[root_i]) 
            {
              candidate_dists_[root_i] = this_weight;
              nearest_neighbors_[root_i] = Edge(i,j,this_weight);
            }
            
          }
         
        } // for j
      } // for i
      
      AddEdges_();
      
    } // while we haven't finished the tree
    
    
  }
  
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::AddEdges_()
  {
    
    // iterate through all the potential edges
    // add them one at a time to avoid cycles
    // DO NOT parallelize naively
    for (int i = 0; i < data_.n_points(); i++)
    {
      
      size_t root_i = components_.Find(i);
      
      Edge& this_edge = nearest_neighbors_[root_i];
      
      size_t root_v = components_.Find(this_edge.v);
      size_t root_u = components_.Find(this_edge.u);
      
      // one of these is always false, but if they're both false, that means
      // the edge is redundant
      if (root_u != root_v) {
        
        // add the edge
        components_.Union(root_u, root_v);
        edge_list_.push_back(this_edge);
        
      }
      
    }
    
    // Don't forget to reset for the next iteration
    for (int i = 0; i < data_.n_points(); i++) {
      
      candidate_dists_[i] = DBL_MAX;
      
    }
    
  } // AddEdges_
  
    
  
  // Single-tree Boruvka, loop over all queries
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::ComputeSTB()
  {
    
    while (edge_list_.size() < data_.n_points() - 1)
    {
    
      // search for each query
      for (size_t i = 0; i < data_.n_points(); i++)
      {
    
        Point<T>& q = data_[i];
        size_t root_q = components_.Find(i);
      
        SearchTree_(q, tree_);
      
      } // loop over queries
    
      AddEdges_();
    
    }
    
  } // ComputeSTB
  
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::SearchTree_(const Point<T>& q,
                                                                  size_t root_q,
                                                                  TTreeType* node)
  {
    
    // we're all connected, so don't search any more
    if (root_q == node->Bound().component()) {
      return;
    }
    else if (EdgePolicy::CanPrune(q, node->Bound(), candidate_dists_[root_q])) {
      return; // we pruned based on distance
    }
    else if (node->IsLeaf())
    {
      for (size_t i = node->Begin(); i < node->End(); i++)
      {
        const Point<T>& point_i = data_[i];
        double this_weight = EdgePolicy::EdgeWeight(q, point_i);
        
        if (this_weight < candidate_dists_[root_q]) 
        {
          
          candidate_dists_[root_q] = this_weight;
          nearest_neighbors_[root_q] = Edge(q, i, this_weight);
          
        }
        
      } // for i
    } // base case
    else {
      
      double left_weight = EdgePolicy::EdgeWeight(q, node->Left()->Center());
      double right_weight = EdgePolicy::EdgeWeight(q, node->Right()->Center());
      
      if (left_weight < right_weight) 
      {
        SearchTree_(q, root_q, node->Left());
        SearchTree_(q, root_q, node->Right());
      }
      else {
        SearchTree_(q, root_q, node->Right());
        SearchTree_(q, root_q, node->Left());
      }
      
      
    } // recursing 
    
  } // SearchTree_()
  
  template<typename T, class EdgePolicy, class TTreeType>
  std::vector<Edge>& MinimumSpanningTree<T, EdgePolicy, TTreeType>::EdgeList() 
  {
    return edge_list_;
  }
  
  
} // namespace



#endif
