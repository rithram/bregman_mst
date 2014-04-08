#include "data.hpp"
#include "KLDivergence.hpp"
#include "L2Divergence.hpp"

using namespace bmst;

int main(int argc, char* argv[])
{
  double eps = 1e-8;
  
  std::vector<double> x_vec(5, 0.1);
  x_vec[0] = 0.6;
  Point<double> x(x_vec);

  std::vector<double> y_vec(5,0.2);  
  Point<double> y(y_vec);

  std::cout << "Testing KL Divergence.\n";
  
  double div_x_y = KLDivergence<double>::BDivergence(x,y);
  double real_div_x_y = 0.3819085009768876;

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
  
  std::vector<double> a_vec(5, 1.0);
  std::vector<double> b_vec(5, 2.0);
  
  a_vec[3] = 0.0;
  
  Point<double> a(a_vec);
  Point<double> b(b_vec);
  
  double div_a_b = KLDivergence<double>::BDivergence(a, b);
  double real_div_a_b = 3.2274112777602184;
  
  assert(fabs(div_a_b - real_div_a_b) < eps);
  
  a[3] = 1.0;
  b[2] = 0.0;
  
  div_a_b = KLDivergence<double>::BDivergence(a,b);
  
  assert(div_a_b == std::numeric_limits<double>::max());
  // assert(div_a_b == -DBL_MAX);
  
  /// NOTE: Removing this test because KL-divergence should not be computed on 
  // negative valued features
  // a[2] = -1.0;
  // div_a_b = KLDivergence<double>::BDivergence(a,b);
  
  // assert(div_a_b == DBL_MAX);
    
  a[2] = 0.0;
  
  div_a_b = KLDivergence<double>::BDivergence(a,b);
  real_div_a_b = 1.2274112777602184;
  
  assert(fabs(real_div_a_b - div_a_b) < eps);
  
  std::cout << "KL Divergence passed.\n";

  std::cout << "Testing L2 Divergence.\n";
  
  div_x_y = L2Divergence<double>::BDivergence(x,y);
  real_div_x_y = 0.1;

  assert(fabs(div_x_y - real_div_x_y) < eps);
  assert(fabs(div_x_y - L2Divergence<double>::BDivergence(y,x)) < eps);

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
