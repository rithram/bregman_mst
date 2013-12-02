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

    bool x_zero = (fabs(x[i]) < std::numeric_limits<double>::epsilon());
    bool y_zero = (fabs(y[i]) < std::numeric_limits<double>::epsilon());
    
    if (!(!y_zero || x_zero)) // y == 0 should imply x == 0, if it doesn't handle specially
    {
      result = x[i] > 0 ? -DBL_MAX : DBL_MAX;
    }
    else if (x_zero) // then they're both zero 
    {
      // 0 log 0 is 0 for KL divergence
      result += y[i];
    }
    else {
      result += x[i] * log(x[i]/y[i]) + y[i] - x[i];      
    }
    
  } // loop over features
    
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
