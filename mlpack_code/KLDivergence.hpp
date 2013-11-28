/**
 *  Implments the KL divergence 
 */

#ifndef KL_DIVERGENCE_HPP_
#define KL_DIVERGENCE_HPP_

#include "data.hpp"
#include <cmath>

namespace bmst {

  /**
   * For now we assume that the distribution is discrete, and represented by a vector of doubles.
   * This means that the user who instantiates this as TBregmanDiv needs to instantiate 
   * TDataType as std::vector<double>
   */
  template<typename T>
  class KLDivergence
  {
    
  public:
    
    static double Divergence(const Point<T>& x, const Point<T>& y);
    
    static Point<T> Gradient(const Point<T>& x);
  
    static Point<T> GradientConjugate(const Point<T>& x);
    
  }; 

}

#include "KLDivergence_impl.hpp"


#endif
