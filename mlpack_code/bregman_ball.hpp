//
//  bregman_ball.hpp
//  
//
//  Created by march on 11/12/13.
//
//

#ifndef _bregman_ball_hpp
#define _bregman_ball_hpp

#include <mlpack/core.hpp>

namespace bbtree {
  
  template <typename TDataType, class TBregmanDiv>
  class BregmanBall <TDataType, TBregmanDiv>
  {
    
  private:
    
    // The center and radius of the ball
    TDataType centroid_;

    double centroid_prime_;
    
    double radius_;
    
  public:
    
    BregmanBall(TDataType& center, double radius);
    
    ~BregmanBall();
    
    // Computes minimum d(x, point) for x in the ball
    double MinRightDist(TDataType& point);
    
    // Computes minimum d(point, x) for x in the ball
    double MinLeftDist(TDataType& point);
    
    // computes minimum d(q,r) with q in this ball and r in the other
    double MinLeftDist(BregmanBall& other);
    
    // computes minimum d(q,r) with r in this ball and q in the other
    double MinRightDist(BregmanBall& other);
    
    double CanPruneRight(TDataType& q, double d_x_c_q);
  
    double CanPruneRight(double theta_l, double theta_r, TDataType& q,
                         double d_x_c_q);
  
    const TDataType& centroid() const;
    
    double centroid_prime() const;
    
    double radius() const;
    
    
  }; // class
  
} // namespace


#include "bregman_ball_impl.hpp"


#endif
