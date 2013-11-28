/**
 *
 */

//
//  bregman_ball_tree_impl.hpp
//  
//
//  Created by march on 11/12/13.
//
//

#ifndef BMST_BREGMAN_BALL_TREE_IMPL_HPP_
#define BMST_BREGMAN_BALL_TREE_IMPL_HPP_

#include "bregman_ball_tree.hpp"

namespace bmst {
  
template <typename T, class TBregmanDiv, class TSplitter>
BregmanBallTree<T, TBregmanDiv, TSplitter>::BregmanBallTree(
    const size_t begin,
    const size_t count,
    const BregmanBall<T, TBregmanDiv>& bounding_ball) :
  begin_(begin),
  count_(count),
  end_(begin + count),
  bounding_ball_(bounding_ball)
{}

template <typename T, class TBregmanDiv, class TSplitter>
BregmanBallTree<T, TBregmanDiv, TSplitter>::BregmanBallTree(
    const Table<T>& table,
    const size_t begin,
    const size_t count) : 
  begin_(begin),
  count_(count),
  end_(begin + count)
{
  // use the table to compute the bounding ball (mean + radius)
  Point<T> center;
  center.zeros(table[begin_].n_dims());
  for (size_t i = begin_; i < end_; i++)
    center += table[i];

  center /= (T) count_;
  double radius = ComputeNodeRadius(table, begin_, end_, center);
  BregmanBall<T, TBregmanDiv> node_bball(center, radius);
  bounding_ball_ = node_bball;
}

template <typename T, class TBregmanDiv, class TSplitter>
void BregmanBallTree<T, TBregmanDiv, TSplitter>::BuildTree(
    Table<T>& data,
    const size_t leaf_size, 
    const double min_ball_width, 
    std::queue<BregmanBallTree<T, TBregmanDiv, TSplitter>*>& node_queue,
    std::vector<size_t>& old_from_new)
{
  typedef BregmanBallTree<T, TBregmanDiv, TSplitter> TNode;
  while (not node_queue.empty()) 
  {
    TNode* current_node = node_queue.front();
    node_queue.pop();

    // std::cout << "Current node count: " << current_node->count_ << 
    //   ", begin @ " << current_node->begin_ << ", end @ " << 
    //   current_node->end_ << std::endl;

    // Try to partition the set:
    // NOTE: Currently, to save one pass over the data, we will always
    // attempt to split the root (and use the left and right stats to 
    // compute the root center)
    std::vector<size_t> membership;
    std::vector<Point<T> > centers;
    std::vector<double> radii;
    TSplitter data_splitter;
    data_splitter.PartitionData(
        data, 
        current_node->begin_,
        current_node->end_, 
        membership,
        centers,
        radii);

    assert(centers.size() == radii.size());
    assert(centers.size() == 2);
    size_t left_count = MatrixSwap(
        data, 
        current_node->begin_, 
        current_node->end_, 
        membership, 
        old_from_new);

    // do something special for the root node
    if (current_node->count_ == data.n_points()) 
    {
      Point<T> root_center 
        = ((double) left_count * centers[0] 
           + (double) (current_node->count_ - left_count) * centers[1]) 
        / (double) current_node->count_;
      double root_radius = 
        ComputeNodeRadius(data, 0, data.n_points(), root_center);
      // initialize the root bounding ball
      BregmanBall<T, TBregmanDiv> root_bball(root_center, root_radius);
      current_node->bounding_ball_ = root_bball;
    }
    if (left_count > 0 and left_count < current_node->Count()) 
    {
      // did find a viable split
      BregmanBall<T, TBregmanDiv> left_bball(centers[0], radii[0]);
      current_node->left_.reset(
          new TNode(current_node->Begin(), left_count, left_bball));
      BregmanBall<T, TBregmanDiv> right_bball(centers[1], radii[1]);
      current_node->right_.reset(new TNode(
          current_node->begin_ + left_count, 
          current_node->count_ - left_count, 
          right_bball));

      // queueing up the children nodes for further tree construction
      if (leaf_size > 0) 
      {
        assert(min_ball_width == 0);
        if (left_count > leaf_size)
          node_queue.push(current_node->left_.get());

        if (current_node->count_ - left_count > leaf_size)
          node_queue.push(current_node->right_.get());
      }
      else 
      {
        assert(min_ball_width > 0);
        assert(leaf_size == 0);
        if (radii[0] > min_ball_width / 2.)
          node_queue.push(current_node->left_.get());

        if (radii[1] > min_ball_width / 2.)
          node_queue.push(current_node->right_.get());
      }
    } // if some split found
  } // node queue loop
} // BuildTree

template <typename T, class TBregmanDiv, class TSplitter>
double BregmanBallTree<T, TBregmanDiv, TSplitter>::ComputeNodeRadius(
    const Table<T>& data,
    const size_t node_begin,
    const size_t node_end,
    const Point<T>& node_center)
{
  double node_radius = 0;
  double div_to_center;
  for (size_t i = node_begin; i < node_end; i++) 
  {
    div_to_center = TBregmanDiv::Divergence(data[i], node_center);
    if (div_to_center > node_radius)
      node_radius = div_to_center;
  }
  return node_radius;
}

template <typename T, class TBregmanDiv, class TSplitter>
size_t BregmanBallTree<T, TBregmanDiv, TSplitter>::MatrixSwap(
    Table<T>& table,
    const size_t node_begin,
    const size_t node_end,
    std::vector<size_t>& membership,
    std::vector<size_t>& old_from_new) 
{
  assert(membership.size() == node_end - node_begin);
  size_t left_ind = 0;
  size_t right_ind = membership.size() - 1;

  while (left_ind < right_ind) 
  {
    while (membership[left_ind] == 0)
      left_ind++;

    while (membership[right_ind] == 1)
      right_ind--;

    if (left_ind > right_ind) 
      break;

    Point<T> temp_point = table[node_begin + left_ind];
    table[node_begin + left_ind] = table[node_begin + right_ind];
    table[node_begin + right_ind] = temp_point;

    size_t temp_ind = old_from_new[node_begin + left_ind];
    old_from_new[node_begin + left_ind] = old_from_new[node_begin + right_ind];
    old_from_new[node_begin + right_ind] = temp_ind;

    assert(membership[left_ind] == 1);
    membership[left_ind] = 0;
    assert(membership[right_ind] == 0);
    membership[right_ind] = 1;
  }
  assert(left_ind == right_ind + 1);
  return left_ind;
}

template <typename T, class TBregmanDiv, class TSplitter>
BregmanBallTree<T, TBregmanDiv, TSplitter>::BregmanBallTree(
    Table<T>& data, 
    std::vector<size_t>& old_from_new,
    const size_t leaf_size, 
    const double min_ball_width) :
  begin_(0),
  count_(data.n_points()),
  end_(data.n_points())
{
  if (leaf_size > 0 and min_ball_width > 0) 
  {
    std::cout << "[ERROR] Can set only one of 'maximum leaf size' and " <<
      "'minimum leaf diameter'" << std::endl;
    exit(1);
  }

  old_from_new.resize(data.n_points());
  for (size_t i = 0; i < count_; i++) 
    old_from_new[i] = i;

  std::queue<BregmanBallTree<T, TBregmanDiv, TSplitter>*> node_queue;
  node_queue.push(this);
  BuildTree(data, leaf_size, min_ball_width, node_queue, old_from_new);
}

template <typename T, class TBregmanDiv, class TSplitter>
BregmanBallTree<T, TBregmanDiv, TSplitter>::~BregmanBallTree()
{
  // nothing to do here since the std::unique_ptr<> should take care 
  // of the automatic deletion of left_ and right_
}

template <typename T, class TBregmanDiv, class TSplitter>
bool BregmanBallTree<T, TBregmanDiv, TSplitter>::CanPruneRight(
    const Point<T>& q,
    const double div_to_center)
{
  // The option where this function can be called with 
  // the divergence(q, centroid) if it is already computed 
  if (div_to_center < std::numeric_limits<double>::max())
    // call bounding_ball CanPruneRight function
    return bounding_ball_.CanPruneRight(q, div_to_center);
  else
    // Compute right divergence to node center
    // call bounding_ball CanPruneRight function
    return bounding_ball_.CanPruneRight(
        q, TBregmanDiv::Divergence(q, bounding_ball_.centroid()));
}

template <typename T, class TBregmanDiv, class TSplitter>
bool BregmanBallTree<T, TBregmanDiv, TSplitter>::CanPruneLeft(
    const Point<T>& q,
    const double div_to_center)
{
  // TO BE FILLED IN (always return false for now)
  return false;
}

template <typename T, class TBregmanDiv, class TSplitter>
bool BregmanBallTree<T, TBregmanDiv, TSplitter>::CanPruneRight(
    const BregmanBallTree<T, TBregmanDiv, TSplitter>& other_node,
    const double node_div_to_center)
{
  // TO BE FILLED IN (always return false for now)
  return false;
}

template <typename T, class TBregmanDiv, class TSplitter>
bool BregmanBallTree<T, TBregmanDiv, TSplitter>::CanPruneLeft(
    const BregmanBallTree<T, TBregmanDiv, TSplitter>& other_node,
    const double center_div_to_node)
{
  // TO BE FILLED IN (always return false for now)
  return false;
}
    
} // namespace

#endif 
