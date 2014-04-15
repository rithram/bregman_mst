/**
 * @file enhanced_bregman_ball_impl.hpp
 *
 * Implementation for EnhancedBregmanBall class
 */

#ifndef BMST_ENHANCED_BREGMAN_BALL_IMPL_HPP_
#define BMST_ENHANCED_BREGMAN_BALL_IMPL_HPP_

#include "enhanced_bregman_ball.hpp"
#include "L2Divergence.hpp"

namespace bmst {

template <typename T, class TBDiv>
EnhancedBregmanBall<T, TBDiv>::EnhancedBregmanBall() :
  TBase(),
  l2_radius_(0),
  jbdiv_radius_(0)
{}

template <typename T, class TBDiv>
EnhancedBregmanBall<T, TBDiv>::EnhancedBregmanBall(
    const Point<T>& right_center, const double right_radius) :
  TBase(right_center, right_radius)
{}

template <typename T, class TBDiv>
EnhancedBregmanBall<T, TBDiv>::EnhancedBregmanBall(
    const Point<T>& right_center, 
    const double right_radius, 
    const Point<T>& left_center, 
    const double left_radius) :
  TBase(right_center, right_radius, left_center, left_radius)
{}

template <typename T, class TBDiv>
EnhancedBregmanBall<T, TBDiv>::~EnhancedBregmanBall()
{}

template <typename T, class TBDiv>
void EnhancedBregmanBall<T, TBDiv>::AddExtraStats(
    const Table<T>& data, const size_t start, const size_t end)
{
  double max_sq_l2_dist = 0;
  double max_sq_jbdiv = 0;
  for (size_t i = 0; i < data.n_points(); ++i) {
    assert(TBase::right_centroid_.n_dims() == data[i].n_dims());
    double sq_l2_dist = L2Divergence<T>::BDivergence(data[i], TBase::right_centroid_);
    if (sq_l2_dist > max_sq_l2_dist)
      max_sq_l2_dist = sq_l2_dist;
    double sq_jbdiv = TBDiv::JBDivergence(data[i], TBase::right_centroid_);
    if (sq_jbdiv > max_sq_jbdiv)
      max_sq_jbdiv = sq_jbdiv;
  }

  l2_radius_ = std::sqrt(max_sq_l2_dist);
  jbdiv_radius_ = std::sqrt(max_sq_jbdiv);

  return;
}

template<typename T, class TBDiv>
bool EnhancedBregmanBall<T, TBDiv>::CanPruneRight(
    const Point<T>& q, const Point<T>& q_prime, const double q_div_to_best_candidate) const
{
  // try pruning using strong convexity
  if (TBDiv::StrongConvexityCoefficient() > 0) 
  {
    const double l2_q_mu = std::sqrt(L2Divergence<T>::BDivergence(q, TBase::right_centroid_));
    if (l2_q_mu > l2_radius_) 
    {
      const double diff = l2_q_mu - l2_radius_;
      const double lb = TBDiv::StrongConvexityCoefficient() * diff * diff;
      if (lb >= q_div_to_best_candidate) 
        return true;
    }
  }

  // check here for CPD-ness, which is known for each divergence
  if (TBDiv::IsCPD()) 
  {
    // try pruning using CPD condition
    const double jbdiv_q_mu = std::sqrt(TBDiv::JBDivergence(q, TBase::right_centroid_));
    if (jbdiv_q_mu > jbdiv_radius_) 
    {
      const double diff = jbdiv_q_mu - jbdiv_radius_;
      const double lb = diff * diff;
      if (lb >= q_div_to_best_candidate)
        return true;
    }
  }

  double d_q_mu = TBDiv::BDivergence(q, TBase::right_centroid_);
  //std::cout << "d_q_mu: " << d_q_mu << ", radius: " << radius_ << 
  // ", q_div_to_best_candidate: " << q_div_to_best_candidate << "\n";
  return TBase::CanPruneRight(q, q_prime, q_div_to_best_candidate, d_q_mu);
}

} // namespace

#endif
