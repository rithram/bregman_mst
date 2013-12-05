/**
 *  Implments the KL divergence 
 */

#ifndef KL_DIVERGENCE_HPP_
#define KL_DIVERGENCE_HPP_

#include "data.hpp"
#include <cmath>

namespace bmst {

template<typename T>
class KLDivergence
{
public:
  static double Divergence(const Point<T>& x, const Point<T>& y);
  static Point<T> Gradient(const Point<T>& x);
  static Point<T> GradientConjugate(const Point<T>& x);
}; // class KLDivergence 

} // namespace

#include "KLDivergence_impl.hpp"


#endif
