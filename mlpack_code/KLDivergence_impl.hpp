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
    if (x[i] < 0 or y[i] < 0) 
    {
      std::cout << "[ERROR] KL divergence cannot be computed for negative "
        "valued features." << std::endl;
      exit(1);
    }

    bool x_zero = (fabs(x[i]) < std::numeric_limits<T>::epsilon());
    bool y_zero = (fabs(y[i]) < std::numeric_limits<T>::epsilon());
    
    if (not x_zero and y_zero) // y == 0 and x > 0,  handle specially
    {
      // x log x - x log y + y - x
      // log 0 = - infty so result = infty
      result = std::numeric_limits<T>::max();
      break; // not need to loop over other features anymore, nothing will reduce infty
    } 
    else if (x_zero and not y_zero) // y > 0, x == 0 , x log (x/y) + y - x = y
    {
      // 0 log 0 is 0 for KL divergence
      result += y[i];
    } 
    else if (not x_zero and not y_zero) 
    {
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
    if (x[i] < 0) 
    {
      std::cout << "[ERROR] Gradient corresponding to KL divergence cannot "
        "be computed for negative valued features." << std::endl;
      exit(1);
    }

    if (x[i] == 0)
      result[i] = -std::numeric_limits<T>::max();
    else
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
