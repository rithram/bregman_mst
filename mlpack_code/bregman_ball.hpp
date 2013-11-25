//
//  bregman_ball.hpp
//  
//
//  Created by march on 11/12/13.
//
//

#ifndef _bregman_ball_hpp
#define _bregman_ball_hpp

#include "data.hpp"

namespace bmst {
  
  template <typename T, class TBregmanDiv>
  class BregmanBall <T, TBregmanDiv>
  {
    
  private:
    
    // The center and radius of the ball
    Point<T> centroid_;

    // The gradient lives in the same space as the data, since we take the 
    // dot product of it with (x - y) in d_f
    Point<T> centroid_prime_;
    
    double radius_;

    bool CanPruneRight(double theta_l, double theta_r, Point<T>& q,
                       double d_x_c_q);
      
  public:
    
    BregmanBall(Point<T>& center, double radius);
    
    ~BregmanBall();
    
    bool CanPruneRight(Point<T>& q, double d_x_c_q);
  
    const Point<T>& centroid() const;
    
    const Point<T>& centroid_prime() const;
    
    double radius() const;
    
    
  }; // class
  
} // namespace


#include "bregman_ball_impl.hpp"


#endif
