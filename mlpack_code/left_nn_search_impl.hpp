
namespace bbtree {

  template<typename TDataType, class TBregmanDiv>
  LeftNNSearch<TDataType, TBregmanDiv>::LeftNNSearch(std::vector<TDataType>& data, int leaf_size, double min_width)
  :
  data_(data),
  leaf_size_(leaf_size),
  min_width_(min_width),
  neighbor_index_(-1),
  neighbor_distance_(DBL_MAX)
  {
  
    tree_ = ConstructBBTree<TDataType, TBregmanDiv>(data, leaf_size_, min_width_, old_from_new_indices_);
    
  }
  
  
  template<typename TDataType, class TBregmanDiv>
  void LeftNNSearch<TDataType, TBregmanDiv>::Compute(TDataType& query) 
  {
    
    neighbor_index_ = -1;
    neighbor_distance_ = DBL_MAX;
    
    double d_x_c_q = TBregmanDiv::Divergence(tree_->centroid(), query);
    
    SearchNode_(tree_, query, d_x_c_q);
    
  }
  
  template<typename TDataType, class TBregmanDiv>
  void LeftNNSearch<TDataType, TBregmanDiv>::SearchNode_(BregmanBallTree<TDataType, TBregmanDiv>* node, TDataType& query, double d_x_c_q) 
  {
    
    if (neighbor_distance_ < DBL_MAX && node->Bound().CanPruneRight(query, d_x_c_q)) 
    {
      
    }
    else if (node->IsLeaf()) {
      
      for (int i = node->Begin(); i < node->End(); i++)
      {
        
        TDataType& ref = data_[i];
        double dist = TBregmanDiv::Divergence(ref, query);
        if (dist < neighbor_distance_) {
          neighbor_distance_ = dist;
          neighbor_index_ = i;
        }
        
      } // for references
      
    } // base case
    else {
      
      double d_left = TBregmanDiv::Divergence(tree_->Left()->centroid(), query);
      double d_right = TBregmanDiv::Divergence(tree_->Right()->centroid(), query);
      
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
