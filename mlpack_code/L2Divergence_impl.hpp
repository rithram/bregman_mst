
namespace bbtree {

  static double L2Divergence<TDataType>::Divergence(TDataType& x, TDataType& y)
  {
    
    // \frac{1}{2} \| x - y \|^2_2
    return 0.5 * (x - y).trans() * (x - y);
    
  }


  // I don't think this should really be a double
  static double L2Divergence<TDataType>::Gradient(TDataType& x)
  {
    // 
    return norm(x);
  }

  static TDataType L2Divergence<TDataType>::GradientConjugate(double x);

  
}