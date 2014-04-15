/**
 * @file enhanced_bregman_ball.hpp
 *
 * This class implements an enhanced version of the bregman ball
 * which contains some extra information (such as L2 radius, 
 * JBDivergence radius, etc) and uses these information and 
 * some assumptions to have new pruning rules
 */

#ifndef BMST_ENHANCED_BREGMAN_BALL_HPP_
#define BMST_ENHANCED_BREGMAN_BALL_HPP_

#include "bregman_ball.hpp"
#include "data.hpp"

namespace bmst {

template <typename T, class TBDiv>
class EnhancedBregmanBall : public BregmanBall<T, TBDiv>
{
private:
  typedef BregmanBall<T, TBDiv> TBase;

  // extra information:
  // max_x sqrt((1/2) ||x - right_centroid_||^2_2)
  double l2_radius_;
  // max_x sqrt(JBDiv(x, right_centroid_))
  double jbdiv_radius_;

public:
  EnhancedBregmanBall();
  EnhancedBregmanBall(const Point<T>& right_center, const double right_radius);
  EnhancedBregmanBall(
      const Point<T>& right_center, 
      const double right_radius, 
      const Point<T>& left_center, 
      const double left_radius);
  
  ~EnhancedBregmanBall();
  
  // Add extra stats from the data if wanted
  // In EnhancedBregmanBall, we compute l2_radius_ and nothing is done here
  void AddExtraStats(const Table<T>& data, const size_t start, const size_t end);

  // Pruning rule for a single query
  bool CanPruneRight(
      const Point<T>& q, 
      const Point<T>& q_prime,
      const double q_div_to_best_candidate) const;
  
}; // class

} // namespace

#endif

#include "enhanced_bregman_ball_impl.hpp"

