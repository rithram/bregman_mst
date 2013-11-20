
#include "KLDivergence.hpp"


using namespace bbtree;


static double KLDivergence::Divergence(std::vector<double>& x, std::vector<double>& y)
{
  
  double result = 0.0;
  
  for (int i = 0; i < x.size(); i++)
  {
    result += x[i] * log(x[i]/y[i]);
  }
  
}

static std::vector<double> KLDivergence::Gradient(std::vector<double>& x)
{

  TDataType result(x.size());

  for (int  i = 0; i < x.size(); i++)
  {
    result[i] = log(x[i]) + 1.0;
  }
  
}

static std::vector<double> KLDivergence::GradientConjugate(std::vector<double>& x)
{
  
  
  
}



