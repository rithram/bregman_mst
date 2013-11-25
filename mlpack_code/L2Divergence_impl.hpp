
#ifndef L2DIVERGENCE_IMPL_HPP_
#define L2DIVERGENCE_IMPL_HPP_

namespace bmst {

  template<typename T>
  double L2Divergence<T>::Divergence(Point<T>& x, Point<T>& y)
  {
    
    // \frac{1}{2} \| x - y \|^2_2
    Point<T> x_minus_y = x - y;
    return 0.5 * Dot(x_minus_y, x_minus_y);
    
  }


  template<typename T>
  Point<T> L2Divergence<T>::Gradient(Point<T>& x)
  {

    return x;

  }

  template<typename T>
  Point<T> L2Divergence<T>::GradientConjugate(Point<T>& x)
  {
    return x;
  }

}

#endif
