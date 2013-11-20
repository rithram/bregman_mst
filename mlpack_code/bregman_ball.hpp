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

    // The gradient lives in the same space as the data, since we take the 
    // dot product of it with (x - y) in d_f
    TDataType centroid_prime_;
    
    double radius_;

    bool CanPruneRight(double theta_l, double theta_r, TDataType& q,
                       double d_x_c_q);
      
  public:
    
    BregmanBall(TDataType& center, double radius);
    
    ~BregmanBall();
    
    bool CanPruneRight(TDataType& q, double d_x_c_q);
  
    const TDataType& centroid() const;
    
    const TDataType& centroid_prime() const;
    
    double radius() const;
    
    
  }; // class
  
} // namespace


#include "bregman_ball_impl.hpp"


#endif
