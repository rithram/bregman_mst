#ifndef BMST_LEFT_NN_SEARCH_HPP_
#define BMST_LEFT_NN_SEARCH_HPP_

#include "bregman_ball_tree.hpp"
#include "kmeans_splitter.hpp"

namespace bmst {

template<typename T, class TBDiv, class TBBall>
class LeftNNSearch {
public:
  
  LeftNNSearch(const Table<T>& data, const size_t leaf_size);
  
  ~LeftNNSearch();
  
  size_t ComputeNeighbor(const Point<T>& query);
  
  size_t ComputeNeighborNaive(const Point<T>& query);
  
private:
  
  typedef KMeansSplitter<T, TBDiv> TSplitter;
  
  typedef BregmanBallTree<T, TBDiv, TBBall, TSplitter> TTreeType;

  Table<T> data_;
  
  TTreeType* tree_;

  size_t leaf_size_;
  
  size_t neighbor_index_;
  double neighbor_distance_;

  std::vector<size_t> old_from_new_indices_;
  
  // functions
  void SearchNode_(
      const TTreeType* node,
      const Point<T>& query,
      const Point<T>& query_prime,
      const T& d_q_to_centroid);
}; // class

} // namespace

#endif

#include "left_nn_search_impl.hpp"


