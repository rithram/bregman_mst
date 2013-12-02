
namespace bmst {

  template<typename T, class TBregmanDiv>
  LeftNNSearch<T, TBregmanDiv>::LeftNNSearch(const Table<T>& data, const size_t leaf_size)
  :
  data_(data),
  leaf_size_(leaf_size),
  neighbor_index_(-1),
  neighbor_distance_(DBL_MAX)
  {
  
    tree_ = new TTreeType(data_, old_from_new_indices_, leaf_size_);
    
  }
  
  template<typename T, class TBregmanDiv>
  LeftNNSearch<T, TBregmanDiv>::~LeftNNSearch()
  {
    
    if (tree_)
      delete tree_;
    
  }
  
  template<typename T, class TBregmanDiv>
  size_t LeftNNSearch<T, TBregmanDiv>::ComputeNeighbor(const Point<T>& query) 
  {
    
    neighbor_index_ = -1;
    neighbor_distance_ = DBL_MAX;
    
    const T& dist_to_centroid = TBregmanDiv::Divergence(query, tree_->Center());
    
    SearchNode_(tree_, query, dist_to_centroid);
    
    return old_from_new_indices_[neighbor_index_];
    
  }
  
  template<typename T, class TBregmanDiv>
  size_t LeftNNSearch<T, TBregmanDiv>::ComputeNeighborNaive(const Point<T>& query)
  {
    
    neighbor_index_ = -1;
    neighbor_distance_ = DBL_MAX;
    
    for (int r = 0; r < data_.n_points(); r++)
    {
      
      double this_dist = TBregmanDiv::Divergence(data_[r], query);
      
      if (this_dist < neighbor_distance_) 
      {
        neighbor_index_ = r;
        neighbor_distance_ = this_dist;
      }
      
    } // loop over references
    
    return old_from_new_indices_[neighbor_index_];
    
  }
  
  template<typename T, class TBregmanDiv>
  void LeftNNSearch<T, TBregmanDiv>::SearchNode_(const TTreeType* node, const Point<T>& query, const T& dist_to_centroid) 
  {
    
    if (neighbor_distance_ < DBL_MAX && node->Bound().CanPruneRight(query, neighbor_distance_)) 
    {
      // then we pruned, so don't do anything
    }
    else if (node->IsLeaf()) {
      
      for (int i = node->Begin(); i < node->End(); i++)
      {
        
        Point<T>& ref = data_[i];
        double dist = TBregmanDiv::Divergence(ref, query);
        if (dist < neighbor_distance_) {
          neighbor_distance_ = dist;
          neighbor_index_ = i;
        }
        
      } // for references
      
    } // base case
    else {
      
      double d_left = TBregmanDiv::Divergence(tree_->Left()->Center(), query);
      double d_right = TBregmanDiv::Divergence(tree_->Right()->Center(), query);
      
      // Prioritize search by distance to centroid
      if (d_left < d_right)
      {
        
        SearchNode_(node->Left(), query, d_left);
        SearchNode_(node->Right(), query, d_right);
        
      }
      else {
        
        SearchNode_(node->Right(), query, d_right);
        SearchNode_(node->Left(), query, d_left);
        
      } 
        
    } // can't prune internal node
    
  } // SearchNode_() 
  
} // namespace
