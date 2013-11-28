/**
 * Implements the L2 Divergence
 */

#ifndef L2DIVERGENCE_HPP_
#define L2DIVERGENCE_HPP_

#include "data.hpp"

namespace bmst {

  // This doesn't actually need templates, since the data type needs to be in R^d for some d
  // Might be templated by the type I'm using for vectors, though? 
  template<typename T>
  class L2Divergence
  {
  
  public:
    
    // static double Divergence(Point<T>& x, Point<T>& y);
    static double Divergence(const Point<T>& x, const Point<T>& y);
    // static double Divergence(Point<T>& x, const Point<T>& y);
    // static double Divergence(const Point<T>& x, Point<T>& y);
  
    // static Point<T> Gradient(Point<T>& x);
    static Point<T> Gradient(const Point<T>& x);
  
    // static Point<T> GradientConjugate(Point<T>& x);
    static Point<T> GradientConjugate(const Point<T>& x);
  
  }; // class

}

#include "L2Divergence_impl.hpp"

#endif
