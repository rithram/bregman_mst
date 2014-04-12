//
//  bregman_ball_impl.hpp
//  
//
//  Created by march on 11/12/13.
//
//

#ifndef BMST_BREGMAN_BALL_IMPL_HPP_
#define BMST_BREGMAN_BALL_IMPL_HPP_

#include "bregman_ball.hpp"

namespace bmst {

template <typename T, class TBregmanDiv>
BregmanBall<T, TBregmanDiv>::BregmanBall() :
  right_centroid_(Point<T>()),
  right_centroid_prime_(Point<T>()),
  left_centroid_(Point<T>()),
  left_centroid_prime_(Point<T>()),
  right_radius_(0),
  left_radius_(0),
  component_(-1)
{}

template <typename T, class TBregmanDiv>
BregmanBall<T, TBregmanDiv>::BregmanBall(
    const Point<T>& right_center, const double right_radius) :
  right_centroid_(right_center),
  right_radius_(right_radius)
{
  right_centroid_prime_  = TBregmanDiv::Gradient(right_center);
}

template <typename T, class TBregmanDiv>
BregmanBall<T, TBregmanDiv>::BregmanBall(
    const Point<T>& right_center, 
    const double right_radius, 
    const Point<T>& left_center, 
    const double left_radius) :
  right_centroid_(right_center),
  left_centroid_(left_center),
  right_radius_(right_radius),
  left_radius_(left_radius)
{
  right_centroid_prime_  = TBregmanDiv::Gradient(right_center);
  left_centroid_prime_  = TBregmanDiv::Gradient(left_center);
}


template <typename T, class TBregmanDiv>
BregmanBall<T, TBregmanDiv>::~BregmanBall()
{}

template<typename T, class TBregmanDiv>
bool BregmanBall<T, TBregmanDiv>::CanPruneRight(
    const Point<T>& q, 
    const Point<T>& q_prime,
    const double q_div_to_best_candidate) const
{
  assert(q.n_dims() == q_prime.n_dims());
  double d_q_mu = TBregmanDiv::BDivergence(q, right_centroid_);
  //std::cout << "d_q_mu: " << d_q_mu << ", radius: " << radius_ << ", q_div_to_best_candidate: " << q_div_to_best_candidate << "\n";
  return CanPruneRight(q, q_prime, q_div_to_best_candidate, d_q_mu);
}

template<typename T, class TBregmanDiv>
bool BregmanBall<T, TBregmanDiv>::CanPruneRight(
    const Point<T>& q,
    const Point<T>& q_prime,
    const double q_div_to_best_candidate, 
    const double q_div_to_centroid) const
{
  // if the query is in the ball, then theta > 1, so don't recurse
  // also, if the query is closer to the centroid than to its candidate, 
  // we can't prune
  // even with theta = 0, so don't bother
  // Also, if the node is a singleton, just do the base case
  if (right_radius_ < std::numeric_limits<T>::epsilon()
      or q_div_to_centroid <= right_radius_ 
      or q_div_to_centroid < q_div_to_best_candidate)
  {
    return false;
  }  

  // initialize at the extreme values of theta
  return CanPruneRight(0.0, 1.0, q, q_prime, q_div_to_best_candidate);
}

template<typename T, class TBregmanDiv>
bool BregmanBall<T, TBregmanDiv>::CanPruneRight(
    const double theta_l,
    const double theta_r,
    const Point<T>& q,
    const Point<T>& q_prime,
    const double q_div_to_best_candidate) const 
{
  if (1.0 - theta_l < std::numeric_limits<T>::epsilon()) 
  {
    // x_theta = mu, but still appears to be outside the ball
    // Possible reasons:
    // * ball radius very small (which we already check for)
    // * d(x_theta, mu) = \infty, which means mu has zero
    //   This is possible when:
    //   ** mu is just a single point with zero (compute d(mu, x) and return),
    //      so effectively no prune
    //   ** all points have same zero (since right_radius_ != \infty), 
    //      so cannot really prune since we dont know the actual d(p, q)
    // NOTE: The above explanation is only valid of KL-divergence
    // but we still cannot prune so that is that
    return false;
  }
  if (theta_r < std::numeric_limits<T>::epsilon()) 
  {
    // x_theta = q, and we are still in the ball
    // * q is almost in the ball, so do not prune
    return false;
  }

  double theta = 0.5 * (theta_l + theta_r);
 
  //std::cout << "Searching, theta = " << theta << "\n";
  
  //std::cout << "q: ";
  //q.print();
  
  Point<T> x_theta_prime = theta * right_centroid_prime_ + (1.0 - theta) * q_prime;
  Point<T> x_theta = TBregmanDiv::GradientConjugate(x_theta_prime);

  //std::cout << "x_theta: ";
  //x_theta.print(); 

  double d_x_theta_mu = TBregmanDiv::BDivergence(x_theta, right_centroid_);
  double d_x_theta_q = TBregmanDiv::BDivergence(x_theta, q);
  
  double L_theta 
      = d_x_theta_q + theta / (1.0 - theta) * (d_x_theta_mu - right_radius_);
  
  //std::cout << "d_x_theta_mu: " << d_x_theta_mu << ", d_x_theta_q: " << d_x_theta_q << ", L_theta: " << L_theta << "\n";
  //std::cout<< "radius: " << right_radius_ << ", q_div_to_best_candidate: " << q_div_to_best_candidate << "\n\n";
      
  if (L_theta > q_div_to_best_candidate)
  {
    return true;
  }
  else if (d_x_theta_mu <= right_radius_ && d_x_theta_q < q_div_to_best_candidate)
  {
    return false;
  }
  else if (d_x_theta_mu > right_radius_)
  {
    // we're outside the ball
    // if (d_x_theta_q < q_div_to_best_candidate) // worth moving inward
    return CanPruneRight(theta, theta_r, q, q_prime, q_div_to_best_candidate);
    // else

  }
  else 
  {
    // we're inside the ball, move outward
    return CanPruneRight(theta_l, theta, q, q_prime, q_div_to_best_candidate);
  }
}

template<typename T, class TBregmanDiv>
bool BregmanBall<T, TBregmanDiv>::CanPruneRight(
    const BregmanBall<T, TBregmanDiv>& other, 
    const double q_div_to_best_candidate, 
    const double q_div_centroids) const
{
  // if the balls intersect, can't prune
  if (q_div_centroids - right_radius_ - other.radius() < 0)
  {
    return false;    
  }
  
  return CanPruneRight(0, 1.0, other, q_div_to_best_candidate, q_div_centroids);
}

template<typename T, class TBregmanDiv>
bool BregmanBall<T, TBregmanDiv>::CanPruneRight(
    const double theta_l,
    const double theta_r,
    const BregmanBall<T, TBregmanDiv>& other, 
    const double q_div_to_best_candidate, 
    const double q_div_centroids) const
{
  return false;
}

} // namespace

#endif
