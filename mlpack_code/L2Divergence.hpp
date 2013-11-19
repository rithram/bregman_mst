/**
 * Implements the L2 Divergence
 */


namespace bbtree {

  // This doesn't actually need templates, since the data type needs to be in R^d for some d
  // Might be templated by the type I'm using for vectors, though? 
  template<typename TDataType>
  class L2Divergence<TDataType> {
  
  private:
    
    
    
  public:
    

    static double Divergence(TDataType& x, TDataType& y);
  
    static double Gradient(TDataType& x);
  
    static TDataType GradientConjugate(double x);
  
  
  
  }; // class

}

#include "L2Divergence_impl.hpp"
