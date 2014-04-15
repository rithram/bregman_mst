/**
 *  Implments the KL divergence 
 */

#ifndef KL_DIVERGENCE_HPP_
#define KL_DIVERGENCE_HPP_

#include "data.hpp"

namespace bmst {

template<typename T>
class KLDivergence
{
public:
  static inline double BDivergence(const Point<T>& x, const Point<T>& y);
  static inline Point<T> Gradient(const Point<T>& x);
  static inline Point<T> GradientConjugate(const Point<T>& x);
  static inline bool IsCPD() { return true; }
  static inline double JBDivergence(const Point<T>& x, const Point<T>& y);
  static inline double StrongConvexityCoefficient() { return 1.0; }
  static size_t bdiv_counter;
  static size_t grad_counter;
  static size_t grad_con_counter;
  static size_t jbdiv_counter;
}; // class KLDivergence 

} // namespace

#endif

#include "KLDivergence_impl.hpp"
