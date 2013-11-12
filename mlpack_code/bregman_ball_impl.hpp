//
//  bregman_ball_impl.hpp
//  
//
//  Created by march on 11/12/13.
//
//

#ifndef _bregman_ball_impl_hpp
#define _bregman_ball_impl_hpp

namespace bbtree {
  
  template <typename TDataType, class TBregmanDiv>
  BregmanBall<TDataType, TBregmanDiv>(TDataType& center, double radius)
  :
  centroid_(center),
  radius_(radius)
  {}
  
  template <typename TDataType, class TBregmanDiv>
  ~BregmanBall<TDataType, TBregmanDiv>()
  {}
  
  
  template<typename TDataType, class TBregmanDiv>
  double BregmanBall<TDataType, TBregmanDiv>::CanPruneRight(TDataType& q,
                                                            TDataType& x_c,
                                                            double d_x_c_q)
  {

    // initialize at the extreme values of theta
    return CanPruneRight(0.0, 1.0, q, x_c, d_x_c_q);
    
  }
  
  
  template<typename TDataType, class TBregmanDiv>
  double BregmanBall<TDataType, TBregmanDiv>::CanPruneRight(double theta_l,
                                                            double theta_r,
                                                            TDataType& q,
                                                            TDataType& x_c,
                                                            double d_x_c_q)
  {
    
    double theta = 0.5 * (theta_l + theta_r);
    double q_prime = TBregmanDiv::Gradient(q);
    
    x_theta = TBregmanDiv::GradientConjugate(theta * centroid_prime_ + (1.0 - theta)*q_prime);

    double d_x_theta_mu = TBregmanDiv::Divergence(x_theta, centroid_);
    double d_x_theta_q = TBregmanDiv::Divergence(x_theta, q);
    
    double L_theta = d_x_theta_q + theta / (1.0 - theta) * (d_x_theta_mu - radius_);
    
    if (L_theta > d_x_c_q)
    {
      return true;
    }
    else if (d_x_theta_mu < radius_ && d_x_theta_q < d_x_c_q)
    {
      
      return false;
      
    }
    else if (d_x_theta_mu > radius_)
    {
      
      // we're outside the ball, so move inward
      return BinarySearch_(theta_l, theta, q, x_c, d_x_c_q);
      
    }
    else {
      
      return BinarySearch_(theta, theta_r, q, x_c, d_x_c_q);
      
    }
    
  }

  // Computes minimum d(x, point) for x in the ball
  template <typename TDataType, class TBregmanDiv>
  bool BregmanBall<TDataType, TBregmanDiv>::CanPruneRight(TDataType& point)
  {
    
    
  }
  
  // Computes minimum d(point, x) for x in the ball
  template <typename TDataType, class TBregmanDiv>
  double BregmanBall<TDataType, TBregmanDiv>::MinLeftDist(TDataType& point);
  
  // computes minimum d(q,r) with q in this ball and r in the other
  template <typename TDataType, class TBregmanDiv>
  double BregmanBall<TDataType, TBregmanDiv>::MinLeftDist(BregmanBall& other);
  
  // computes minimum d(q,r) with r in this ball and q in the other
  template <typename TDataType, class TBregmanDiv>
  double BregmanBall<TDataType, TBregmanDiv>::MinRightDist(BregmanBall& other);

  
  
}





#endif
