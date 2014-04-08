#ifndef L2DIVERGENCE_IMPL_HPP_
#define L2DIVERGENCE_IMPL_HPP_

#include "L2Divergence.hpp"

namespace bmst {

template<typename T> 
size_t L2Divergence<T>::bdiv_counter = 0;

template<typename T> 
size_t L2Divergence<T>::grad_counter = 0;

template<typename T> 
size_t L2Divergence<T>::grad_con_counter = 0;

template<typename T> 
size_t L2Divergence<T>::jbdiv_counter = 0;

template<typename T>
double L2Divergence<T>::BDivergence(const Point<T>& x, const Point<T>& y)
{
  ++bdiv_counter;  
  // \frac{1}{2} \| x - y \|^2_2
  Point<T> x_minus_y = x - y;
  return 0.5 * Dot(x_minus_y, x_minus_y);
}

template<typename T>
Point<T> L2Divergence<T>::Gradient(const Point<T>& x)
{
  ++grad_counter;
  return x;
}

template<typename T>
Point<T> L2Divergence<T>::GradientConjugate(const Point<T>& x)
{
  ++grad_con_counter;
  return x;
}

template<typename T>
double L2Divergence<T>::JBDivergence(const Point<T>& x, const Point<T>& y)
{
  // compute f(x) + f(x) - 2 f((x+y)/2)
  // 0.5 ||x||^2 + 0.5 ||y||^2 - ||(x + y) / 2||^2
  //  = 0.25 * || x - y ||^2 
  ++jbdiv_counter;
  Point<T> x_minus_y = x - y;
  return 0.25 * Dot(x_minus_y, x_minus_y);
}

}

#endif
