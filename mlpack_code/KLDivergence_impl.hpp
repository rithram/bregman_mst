#ifndef KL_DIVERGENCE_IMPL_HPP_
#define KL_DIVERGENCE_IMPL_HPP_


namespace bmst {

  template<typename T>
  double KLDivergence<T>::Divergence(const Point<T>& x, const Point<T>& y)
  {
  
    assert(x.n_dims() == y.n_dims());
  
    double result = 0.0;
  
    for (int i = 0; i < x.n_dims(); i++)
    {
      result += x[i] * log(x[i]/y[i]) + y[i] - x[i];
    }
    
    return result;
  
  }

  template<typename T>
  Point<T> KLDivergence<T>::Gradient(const Point<T>& x)
  {

    Point<T> result;
    result.zeros(x.n_dims());

    for (int  i = 0; i < x.n_dims(); i++)
    {
      result[i] = log(x[i]) + 1.0;
    }
  
    return result;
  
  }

  template<typename T>
  Point<T> KLDivergence<T>::GradientConjugate(const Point<T>& x)
  {
  
    Point<T> result;
    result.zeros(x.n_dims());
    
    for (int i = 0; i < x.n_dims(); i++)
    {
      result[i] = exp(x[i] - 1.0);
    }
    
    return result;
  
  }


}

#endif
