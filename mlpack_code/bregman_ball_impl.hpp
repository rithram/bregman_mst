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
  {
    
    centroid_prime_  = TBregmanDiv::Gradient(centroid_);
    
  }
  
  template <typename TDataType, class TBregmanDiv>
  ~BregmanBall<TDataType, TBregmanDiv>()
  {}
  
  
  template<typename TDataType, class TBregmanDiv>
  bool BregmanBall<TDataType, TBregmanDiv>::CanPruneRight(TDataType& q,
                                                          double d_x_c_q)
  {

    TDataType q_prime = TBregmanDiv::Gradient(q);

    // initialize at the extreme values of theta
    return CanPruneRight(0.0, 1.0, q, q_prime, d_x_c_q);
    
  }
  
  
  template<typename TDataType, class TBregmanDiv>
  bool BregmanBall<TDataType, TBregmanDiv>::CanPruneRight(double theta_l,
                                                            double theta_r,
                                                            TDataType& q,
                                                            TDataType& q_prime,
                                                            double d_x_c_q)
  {
    
    double theta = 0.5 * (theta_l + theta_r);
    
    TData_type x_theta = TBregmanDiv::GradientConjugate(theta * centroid_prime_ + (1.0 - theta)*q_prime);

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
      return CanPruneRight(theta_l, theta, q, q_prime, d_x_c_q);
      
    }
    else {
      
      // we're inside the ball, move outward
      return CanPruneRight(theta, theta_r, q, q_prime, d_x_c_q);
      
    }
    
  }
  
  template<typename TDataType, class TBregmanDiv>
  const TDataType& BregmanBall<TDataType, TBregmanDiv>::centroid() const
  {
    return centroid_;
  }
  
  template<typename TDataType, class TBregmanDiv>
  const TDataType& BregmanBall<TDataType, TBregmanDiv>::centroid_prime() const
  {
    return centroid_prime_;
  }
  
  template<typename TDataType, class TBregmanDiv>
  double BregmanBall<TDataType, TBregmanDiv>::radius() const
  {
    return radius_;
  }
  
  
}





#endif
