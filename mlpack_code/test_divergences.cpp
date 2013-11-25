
#include "data.hpp"
#include "KLDivergence.hpp"
#include "L2Divergence.hpp"

using namespace bmst;

int main(int argc, char* argv[])
{
  
  double eps = 1e-8;
  
  std::vector<double> x_vec(5, 0.1);
  x_vec[0] = 0.5;
  Point<double> x(x_vec);

  std::vector<double> y_vec(5,0.2);  
  Point<double> y(y_vec);

  std::cout << "Testing KL Divergence.\n";
  
  double div_x_y = KLDivergence<double>::Divergence(x,y);
  double real_div_x_y = 0.18088649371309945;

  assert(fabs(div_x_y - real_div_x_y) < eps);

  Point<double> x_prime = KLDivergence<double>::Gradient(x);
  Point<double> x_prime_inv = KLDivergence<double>::GradientConjugate(x_prime);

  Point<double> y_prime = KLDivergence<double>::Gradient(y);
  Point<double> y_prime_inv = KLDivergence<double>::GradientConjugate(y_prime);
  
  for (int i = 0; i < x.n_dims(); i++)
  {
    assert (fabs(x[i] - x_prime_inv[i]) < eps);
    assert (fabs(y[i] - y_prime_inv[i]) < eps);
  }
  
  std::cout << "KL Divergence passed.\n";

  std::cout << "Testing L2 Divergence.\n";
  
  div_x_y = L2Divergence<double>::Divergence(x,y);
  real_div_x_y = 0.5 * 0.13;

  assert(fabs(div_x_y - real_div_x_y) < eps);
  assert(fabs(div_x_y - L2Divergence<double>::Divergence(y,x)) < eps);

  x_prime = L2Divergence<double>::Gradient(x);
  x_prime_inv = L2Divergence<double>::GradientConjugate(x_prime);

  y_prime = L2Divergence<double>::Gradient(y);
  y_prime_inv = L2Divergence<double>::GradientConjugate(y_prime);
  
  for (int i = 0; i < x.n_dims(); i++)
  {
    assert (fabs(x[i] - x_prime_inv[i]) < eps);
    assert (fabs(y[i] - y_prime_inv[i]) < eps);
  }
  
  std::cout << "L2 Divergence passed.\n";
  




  return 0;
}

