/**
 * Implements the L2 Divergence
 */

#ifndef L2DIVERGENCE_HPP_
#define L2DIVERGENCE_HPP_

#include "data.hpp"

namespace bmst {

template<typename T>
class L2Divergence
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
}; // class L2Divergence

} // namespace

#endif

#include "L2Divergence_impl.hpp"
