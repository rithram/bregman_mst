#ifndef BMST_LEFT_NN_SEARCH_IMPL_HPP_
#define BMST_LEFT_NN_SEARCH_IMPL_HPP_

#include "left_nn_search.hpp"

namespace bmst {

template<typename T, class TBDiv, class TBBall>
LeftNNSearch<T, TBDiv, TBBall>::LeftNNSearch(
    const Table<T>& data, const size_t leaf_size) :
  data_(data),
  leaf_size_(leaf_size),
  neighbor_index_(-1),
  neighbor_distance_(std::numeric_limits<T>::max())
{
  tree_ = new TTreeType(data_, old_from_new_indices_, leaf_size_);
}

template<typename T, class TBDiv, class TBBall>
LeftNNSearch<T, TBDiv, TBBall>::~LeftNNSearch()
{
  if (tree_)
    delete tree_;
}

template<typename T, class TBDiv, class TBBall>
size_t LeftNNSearch<T, TBDiv, TBBall>::ComputeNeighbor(const Point<T>& query) 
{
  neighbor_index_ = -1;
  neighbor_distance_ = std::numeric_limits<T>::max();
  
  const T& dist_to_centroid = TBDiv::BDivergence(query, tree_->RCenter());
  
  SearchNode_(tree_, query, dist_to_centroid);
  
  if (neighbor_index_ == -1) {
    assert(neighbor_distance_ == std::numeric_limits<T>::max());
    return -1;
  } else {
    assert(neighbor_distance_ < std::numeric_limits<T>::max());
    assert(old_from_new_indices_[neighbor_index_]  < data_.n_points());
    return old_from_new_indices_[neighbor_index_];
  }
}

template<typename T, class TBDiv, class TBBall>
size_t LeftNNSearch<T, TBDiv, TBBall>::ComputeNeighborNaive(const Point<T>& query)
{
  neighbor_index_ = -1;
  neighbor_distance_ = std::numeric_limits<T>::max();
  
  for (int r = 0; r < data_.n_points(); r++)
  {
    double this_dist = TBDiv::BDivergence(data_[r], query);
    
    if (this_dist < neighbor_distance_) 
    {
      neighbor_index_ = r;
      neighbor_distance_ = this_dist;
    }
  } // loop over references

  if (neighbor_index_ == -1) {
    assert(neighbor_distance_ == std::numeric_limits<T>::max());
    return -1;
  } else {
    assert(neighbor_distance_ < std::numeric_limits<T>::max());
    assert(old_from_new_indices_[neighbor_index_]  < data_.n_points());
    return old_from_new_indices_[neighbor_index_];
  }
}

template<typename T, class TBDiv, class TBBall>
void LeftNNSearch<T, TBDiv, TBBall>::SearchNode_(
    const TTreeType* node, const Point<T>& query, const T& dist_to_centroid) 
{
  if (neighbor_distance_ < std::numeric_limits<T>::max() 
      and node->Bound().CanPruneRight(query, neighbor_distance_)) 
  {
    // then we pruned, so don't do anything
    return;
  }
  else if (node->IsLeaf()) 
  {
    for (int i = node->Begin(); i < node->End(); i++)
    {
      double dist = TBDiv::BDivergence(data_[i], query);
      if (dist < neighbor_distance_) {
        neighbor_distance_ = dist;
        neighbor_index_ = i;
      }
    } // for references
  } // base case
  else 
  {
    double d_left = TBDiv::BDivergence(tree_->Left()->RCenter(), query);
    double d_right = TBDiv::BDivergence(tree_->Right()->RCenter(), query);
    
    // Prioritize search by distance to centroid
    if (d_left < d_right)
    {
      SearchNode_(node->Left(), query, d_left);
      SearchNode_(node->Right(), query, d_right);
    }
    else 
    {
      SearchNode_(node->Right(), query, d_right);
      SearchNode_(node->Left(), query, d_left);
    } 
  } // can't prune internal node
  return;
} // SearchNode_() 

} // namespace

#endif
