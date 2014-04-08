/**
 * @file 
 */

//
//  bregman_ball_tree.h
//  
//
//  Created by march on 11/12/13.
//
//

#ifndef BMST_BREGMAN_BALL_TREE_HPP_
#define BMST_BREGMAN_BALL_TREE_HPP_

#include <memory>
#include <queue>
#include <vector>

#include "data.hpp"

namespace bmst {

template <typename T, class TBDiv, class TBBall, class TSplitter>
class BregmanBallTree
{
private:
  typedef BregmanBallTree<T, TBDiv, TBBall, TSplitter> TBBTree;

  // Indices into the data set
  int begin_;
  int end_;
  int count_;

  // Node bounds
  TBBall bounding_ball_;

  // children
  std::unique_ptr<TBBTree> left_;
  std::unique_ptr<TBBTree> right_;

  // Initializer
  BregmanBallTree(
      const size_t begin,
      const size_t count,
      const TBBall& bounding_ball);

  // Initializer
  BregmanBallTree(
      const Table<T>& table,
      const size_t begin,
      const size_t count);

  // Helper functions
  void BuildTree(
      Table<T>& table,
      const size_t leaf_size,
      const double min_ball_width,
      std::queue<TBBTree*>& node_queue,
      std::vector<size_t>& old_from_new);

  double ComputeNodeRadius(
      const Table<T>& data,
      const size_t node_begin,
      const size_t node_end,
      const Point<T>& node_center);

  size_t MatrixSwap(
      Table<T>& table,
      const size_t node_begin,
      const size_t node_end,
      std::vector<size_t>& membership,
      std::vector<size_t>& old_from_new);

public:
  // Initializer
  BregmanBallTree(
      Table<T>& data, 
      std::vector<size_t>& old_from_new,
      const size_t leaf_size = 10, 
      const double min_ball_width = 0);
    
  ~BregmanBallTree();
  // Tree info accessors
  bool IsLeaf() const { if (left_) return false; else return true; };
  TBBTree* Left() const { return left_.get(); }
  TBBTree* Right() const { return right_.get(); }
  const int Begin() const { return begin_; }
  const int End() const { return end_; } 
  const int Count() const { return count_; }
    
  // node statistic accessors
  const Point<T>& RCenter() const { return bounding_ball_.right_centroid(); }
  const double RRadius() const { return bounding_ball_.right_radius(); }
  const Point<T>& LCenter() const { return bounding_ball_.left_centroid(); }
  const double LRadius() const { return bounding_ball_.left_radius(); }
  const TBBall& Bound() const { return bounding_ball_; }
  TBBall& Bound() { return bounding_ball_; }

  // Pruning functions
  // point-ball right-prune
  bool CanPruneRight(
      const Point<T>& q, 
      const double q_div_to_best_candidate,
      const double div_to_center = std::numeric_limits<double>::max());

  // TO-DO: if needed, implement point-ball left-prune
  // API might change
  bool CanPruneLeft(
      const Point<T>& q, 
      const double q_div_to_best_candidate,
      const double div_to_center = std::numeric_limits<double>::max());

  // TO-DO: ball-ball right prune
  bool CanPruneRight(
      const TBBTree& other_node,
      const double node_max_div_to_best_candidate,
      const double node_div_to_center = std::numeric_limits<double>::max());

  // TO-DO: if needed, ball-ball left prune
  // API might change
  bool CanPruneLeft(
      const TBBTree& other_node,
      const double node_max_div_to_best_candidate,
      const double center_div_to_node = std::numeric_limits<double>::max());
}; // class
  
} // namespace

#endif

#include "bregman_ball_tree_impl.hpp"

