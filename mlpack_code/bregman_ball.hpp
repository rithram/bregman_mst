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
  Point<T> centroid_;

  // The gradient lives in the same space as the data, since we take the 
  // dot product of it with (x - y) in d_f
  Point<T> centroid_prime_;
  
  double radius_;

  bool CanPruneRight(
      const double theta_l, 
      const double theta_r, 
      const Point<T>& q, 
      const Point<T>& q_prime,
      const double q_div_to_best_candidate) const;
    
public:
  BregmanBall();
  BregmanBall(const Point<T>& center, const double radius);
  
  ~BregmanBall();
  
  bool CanPruneRight(const Point<T>& q, const double q_div_to_best_candidate) const;
  
  // We'll precompute this divergence to prioritize the tree search, so this function allows us not to compute
  // the distance to the centroid again
  bool CanPruneRight(const Point<T>& q, const double q_div_to_best_candidate, const double q_div_to_centroid) const;

  const Point<T>& centroid() const;
  
  const Point<T>& centroid_prime() const;
  
  const double radius() const;
}; // class

} // namespace

#include "bregman_ball_impl.hpp"

#endif
