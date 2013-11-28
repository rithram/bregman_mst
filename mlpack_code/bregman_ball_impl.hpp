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
  centroid_(Point<T>()),
  centroid_prime_(Point<T>()),
  radius_(0)
{}

template <typename T, class TBregmanDiv>
BregmanBall<T, TBregmanDiv>::BregmanBall(Point<T>& center, double radius) :
  centroid_(center),
  radius_(radius)
{
  centroid_prime_  = TBregmanDiv::Gradient(centroid_);
}

template <typename T, class TBregmanDiv>
BregmanBall<T, TBregmanDiv>::~BregmanBall()
{}

template<typename T, class TBregmanDiv>
bool BregmanBall<T, TBregmanDiv>::CanPruneRight(
    Point<T>& q, double d_x_c_q)
{
  // if the query is in the ball, then theta > 1, so don't recurse
  // also, if the query is closer to the centroid than to its candidate, we can't prune
  // even with theta = 0, so don't bother
  double d_q_mu = TBregmanDiv::Divergence(q, centroid_);
  //std::cout << "d_q_mu: " << d_q_mu << ", radius: " << radius_ << ", d_x_c_q: " << d_x_c_q << "\n";
  if (d_q_mu <= radius_ || d_q_mu < d_x_c_q)
  {
    return false;
  }  

  Point<T> q_prime = TBregmanDiv::Gradient(q);

  // initialize at the extreme values of theta
  return CanPruneRight(0.0, 1.0, q, q_prime, d_x_c_q);
}


template<typename T, class TBregmanDiv>
bool BregmanBall<T, TBregmanDiv>::CanPruneRight(
    double theta_l,
    double theta_r,
    Point<T>& q,
    Point<T>& q_prime,
    double d_x_c_q)
{
  double theta = 0.5 * (theta_l + theta_r);
  
  //std::cout << "Searching, theta = " << theta << "\n";
  
  //std::cout << "q: ";
  //q.print();
  
  Point<T> x_theta_prime = theta * centroid_prime_ + (1.0 - theta) * q_prime;
  Point<T> x_theta = TBregmanDiv::GradientConjugate(x_theta_prime);

  //std::cout << "x_theta: ";
  //x_theta.print(); 

  double d_x_theta_mu = TBregmanDiv::Divergence(x_theta, centroid_);
  double d_x_theta_q = TBregmanDiv::Divergence(x_theta, q);
  
  double L_theta = d_x_theta_q + theta / (1.0 - theta) * (d_x_theta_mu - radius_);
  
  //std::cout << "d_x_theta_mu: " << d_x_theta_mu << ", d_x_theta_q: " << d_x_theta_q << ", L_theta: " << L_theta << "\n";
  //std::cout<< "radius: " << radius_ << ", d_x_c_q: " << d_x_c_q << "\n\n";
      
  if (L_theta > d_x_c_q)
  {
    return true;
  }
  else if (d_x_theta_mu <= radius_ && d_x_theta_q < d_x_c_q)
  {
    return false;
  }
  else if (d_x_theta_mu > radius_)
  {
    // we're outside the ball, so move inward
    return CanPruneRight(theta, theta_r, q, q_prime, d_x_c_q);
  }
  else 
  {
    // we're inside the ball, move outward
    return CanPruneRight(theta_l, theta, q, q_prime, d_x_c_q);
  }
}

template<typename T, class TBregmanDiv>
const Point<T>& BregmanBall<T, TBregmanDiv>::centroid() const
{
  return centroid_;
}

template<typename T, class TBregmanDiv>
const Point<T>& BregmanBall<T, TBregmanDiv>::centroid_prime() const
{
  return centroid_prime_;
}

template<typename T, class TBregmanDiv>
double BregmanBall<T, TBregmanDiv>::radius() const
{
  return radius_;
}

} // namespace

#endif
