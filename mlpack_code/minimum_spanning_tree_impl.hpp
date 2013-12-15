
#ifndef MINIMUM_SPANNING_TREE_IMPL_HPP_
#define MINIMUM_SPANNING_TREE_IMPL_HPP_

namespace bmst {

  template<typename T, class EdgePolicy, class TTreeType>
  MinimumSpanningTree<T, EdgePolicy, TTreeType>::MinimumSpanningTree(Table<T>& data, int leaf_size)
  :
  data_(data),
  components_(data.n_points()),
  nearest_neighbors_(data.n_points()),
  candidate_dists_(data.n_points(), DBL_MAX)
  {
    
    tree_ = new TTreeType(data_, old_from_new_, leaf_size);
    
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
    
    // make sure we reset everything first
    ResetAll_();
    
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
          
          if (root_q != root_r and this_weight < candidate_dists_[root_q]) 
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
    
    // make sure we reset everything first
    components_.Reset();
    
    std::vector<std::vector<double> > edge_weights;
    
    if (use_n2_memory)
    {

      edge_weights.resize(data_.n_points());  
      
      // allocate the memory first
      for (size_t i = 0; i < data_.n_points(); i++)
      {
        
        edge_weights[i].resize(data_.n_points());
        
      } // for i
      
      for (size_t j = 0; j < data_.n_points(); j++)
      {
        
        edge_weights[j][j] = 0.0;
        
        for (size_t k = j+1; k < data_.n_points(); k++)
        {
      
          // this assumes that the edge weights are symmetric
          double this_weight = EdgePolicy::EdgeWeight(data_[j], data_[k]);
          
          edge_weights[j][k] = this_weight;
          edge_weights[k][j] = this_weight;
          
        } // for k
      } // for j
      
    } // if we're precomputing all of the distances
    
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
        size_t root_i = components_.Find(i);
        
        for (size_t j = 0; j < data_.n_points(); j++)
        {
        
          // don't bother if they're already connected
          if (i == j || root_i == components_.Find(j)) continue;
        
          const Point<T>& point_j = data_[j];
          
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
      
      if (root_u != root_v) {
        
        // add the edge
        components_.Union(root_u, root_v);
        
        edge_list_.push_back(this_edge);
        
        //std::cout << "Adding edges: (" << this_edge.u << ", " << this_edge.v << "): " << this_edge.weight << "\n\n\n";
        
      }
      
    }
    
    // Don't forget to reset for the next iteration
    for (int i = 0; i < data_.n_points(); i++) {
      
      candidate_dists_[i] = DBL_MAX;
      nearest_neighbors_[i].weight = DBL_MAX;
      
    }
    
  } // AddEdges_
  
    
  
  // Single-tree Boruvka, loop over all queries
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::ComputeSTB()
  {
    
    // Reset in case we already used this object
    ResetAll_();
    
    while (edge_list_.size() < data_.n_points() - 1)
    {
    
      // search for each query
      for (size_t i = 0; i < data_.n_points(); i++)
      {
    
        Point<T>& q = data_[i];
        size_t root_q = components_.Find(i);
      
        SearchTree_(q, i, root_q, tree_);
      
      } // loop over queries
    
      AddEdges_();
      
      UpdateTree_(tree_);
    
    }
    
  } // ComputeSTB
  
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::UpdateTree_(TTreeType* node) 
  {
    
    if (node->IsLeaf()) 
    {
      
      if (node->Bound().Component() < 0) 
      {
        
        size_t comp = components_.Find(node->Begin());
        
        for (size_t i = node->Begin() + 1; i < node->End(); i++) {
          
          size_t this_comp = components_.Find(i);
          
          if (this_comp != comp) {
            comp = -1;
            break;
          }
          
        } // loop over points in the leaf
        
        node->Bound().SetComponent(comp);
        
      } // do we need to see if the node is connected now?
      else {
        
        // it might have changed id, so just set it to the new one
        // There is no way for the node to have stopped being connected
        node->Bound().SetComponent(components_.Find(node->Begin()));
        
      }
      
    } // is the node a leaf?
    else {
      
      UpdateTree_(node->Left());
      UpdateTree_(node->Right());
      
      size_t left_comp = node->Left()->Bound().Component();
      size_t right_comp = node->Right()->Bound().Component();

      // We don't need to check if they're positive, since this is ok in the
      // case that they're both -1      
      if (left_comp == right_comp) 
      {
        node->Bound().SetComponent(left_comp);
      }
      // we assume that it's already -1
      
    }
    
  }
  
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::ResetTree_(TTreeType* node) 
  {
    
    if (not node->IsLeaf())
    {
      ResetTree_(node->Left());
      ResetTree_(node->Right());
    }
    
    node->Bound().SetComponent(-1);
    
  }
  
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::SearchTree_(const Point<T>& q,
                                                                  size_t q_index,
                                                                  size_t root_q,
                                                                  TTreeType* node)
  {
    
    // we're all connected, so don't search any more
    if (root_q == node->Bound().Component()) {
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
          nearest_neighbors_[root_q] = Edge(q_index, i, this_weight);
          
        }
        
      } // for i
    } // base case
    else {
      
      double left_weight = EdgePolicy::EdgeWeight(q, node->Left()->RCenter());
      double right_weight = EdgePolicy::EdgeWeight(q, node->Right()->RCenter());
      
      if (left_weight < right_weight) 
      {
        SearchTree_(q, q_index, root_q, node->Left());
        SearchTree_(q, q_index, root_q, node->Right());
      }
      else {
        SearchTree_(q, q_index, root_q, node->Right());
        SearchTree_(q, q_index, root_q, node->Left());
      }
      
      
    } // recursing 
    
  } // SearchTree_()
  
  template<typename T, class EdgePolicy, class TTreeType>
  std::vector<Edge>& MinimumSpanningTree<T, EdgePolicy, TTreeType>::EdgeList() 
  {
    
    std::sort(edge_list_.begin(), edge_list_.end(), EdgeSorter);
    
    for (int i = 0; i < edge_list_.size(); i++)
    {
      edge_list_[i].u = old_from_new_[edge_list_[i].u];
      edge_list_[i].v = old_from_new_[edge_list_[i].v];
    }
    
    return edge_list_;
  
  }
  
  template<typename T, class EdgePolicy, class TTreeType>
  void MinimumSpanningTree<T, EdgePolicy, TTreeType>::ResetAll_() 
  {
    
    components_.Reset();
    
    // don't call UpdateTree_ here, won't work because leaves already have 
    // their components set
    ResetTree_(tree_);
    
  }
  
  
} // namespace



#endif
