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
  
  const T dist_to_centroid = TBDiv::BDivergence(query, tree_->RCenter());
  const Point<T> query_prime = TBDiv::Gradient(query);
  
  SearchNode_(tree_, query, query_prime, dist_to_centroid);
  
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
    const TTreeType* node, 
    const Point<T>& query, 
    const Point<T>& query_prime, 
    const T dist_to_centroid) 
{
  // at leaf, do exhaustive search
  if (node->IsLeaf()) 
  {
    for (int i = node->Begin(); i < node->End(); i++)
    {
      double dist = TBDiv::BDivergence(data_[i], query);
      if (dist < neighbor_distance_) 
      {
        neighbor_distance_ = dist;
        neighbor_index_ = i;
      }
    } // for references
    return;
  } // base case

  double d_left = TBDiv::BDivergence(tree_->Left()->RCenter(), query);
  double d_right = TBDiv::BDivergence(tree_->Right()->RCenter(), query);
  // Prioritize search by distance to centroid
  // NOTE: This current scheme always goes to atleast one leaf of any subtree 
  // that is not pruned -- this is useful if you are not pruning a lot anyways
  // because in that case, you save computation that is needed for doing the 
  // pruning check
  if (d_left < d_right)
  {
    // search left
    // if (not node->Left()->Bound().CanPruneRight(
    //     query, query_prime, neighbor_distance_))
    SearchNode_(node->Left(), query, query_prime, d_left);
    // try to prune right
    if (not node->Right()->Bound().CanPruneRight(
        query, query_prime, neighbor_distance_))
      SearchNode_(node->Right(), query, query_prime, d_right);
  }
  else 
  {
    // search right
    // if (not node->Right()->Bound().CanPruneRight(
    //     query, query_prime, neighbor_distance_))
    SearchNode_(node->Right(), query, query_prime, d_right);
    // try to prune left
    if (not node->Left()->Bound().CanPruneRight(
        query, query_prime, neighbor_distance_))
      SearchNode_(node->Left(), query, query_prime, d_left);
  } 

  return;
} // SearchNode_() 

} // namespace

#endif
