//
//  bregman_ball.hpp
//  
//
//  Created by march on 11/12/13.
//
//

#ifndef BMST_BREGMAN_BALL_HPP_
#define BMST_BREGMAN_BALL_HPP_

#include "data.hpp"

namespace bmst {

template <typename T, class TBregmanDiv>
class BregmanBall 
{
private:
  // The center and radius of the ball

  // the right centroid is the point such that 
  // Div(x, right_centroid_) \leq right_radius_ for all x in the ball
  Point<T> right_centroid_;

  // the left centroid is the point such that
  // Div(left_centroid_, x) \leq left_radius_ for all x in the ball
  Point<T> left_centroid_;

  // the gradients of the centroids
  Point<T> right_centroid_prime_;
  Point<T> left_centroid_prime_;
  
  // max_x BDiv(x, right_centroid_)
  double right_radius_;
  // max_x BDiv(left_centroid_, x)
  double left_radius_;
  
  size_t component_;

  // helper for pruning in single tree traversal
  bool CanPruneRight(
      const double theta_l, 
      const double theta_r, 
      const Point<T>& q, 
      const Point<T>& q_prime,
      const double q_div_to_best_candidate) const;
    
  bool CanPruneRight(
      const double theta_l,
      const double theta_r,
      const BregmanBall<T, TBregmanDiv>& other, 
      const double q_div_to_best_candidate, 
      const double q_div_centroids) const;

public:
  BregmanBall();
  BregmanBall(const Point<T>& right_center, const double right_radius);
  BregmanBall(
      const Point<T>& right_center, 
      const double right_radius, 
      const Point<T>& left_center, 
      const double left_radius);
  
  ~BregmanBall();

  // Add extra stats from the data if wanted
  // In plain BregmanBall, nothing is done here
  void AddExtraStats(const Table<T>& data, const size_t start, const size_t end) {}
  
  // Pruning rule for a single query
  bool CanPruneRight(
      const Point<T>& q,
      const Point<T>& q_prime,
      const double q_div_to_best_candidate) const;
  
  // We'll precompute this divergence to prioritize the tree search, 
  // so this function allows us not to compute
  // the distance to the centroid again
  bool CanPruneRight(
      const Point<T>& q,
      const Point<T>& q_prime,
      const double q_div_to_best_candidate, 
      const double q_div_to_centroid) const;

  // pruning rule for two nodes, with the query on the left
  bool CanPruneRight(
      const BregmanBall<T, TBregmanDiv>& other, 
      const double q_div_to_best_candidate, 
      const double q_div_centroids) const;

  const Point<T>& left_centroid() const { return left_centroid_; }
  
  const Point<T>& left_centroid_prime() const { return left_centroid_prime_; }

  const Point<T>& right_centroid() const { return right_centroid_; }
  
  const Point<T>& right_centroid_prime() const { return right_centroid_prime_; }
  
  const double right_radius() const { return right_radius_; }

  const double left_radius() const { return left_radius_; }
  
  const size_t Component() const { return component_; }
  
  void SetComponent(size_t comp) { component_ = comp; }
  
}; // class

} // namespace

#endif

#include "bregman_ball_impl.hpp"

