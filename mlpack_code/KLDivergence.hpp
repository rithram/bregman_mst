/**
 *  Implments the KL divergence 
 */

namespace bbtree {

  /**
   * For now we assume that the distribution is discrete, and represented by a vector of doubles.
   * This means that the user who instantiates this as TBregmanDiv needs to instantiate 
   * TDataType as std::vector<double>
   */
  class KLDivergence
  {
    
  private:
    
    
    
    
  public:
    
    static double Divergence(std::vector<double>& x, std::vector<double>& y);
    
    static std::vector<double> Gradient(std::vector<double>& x);
  
    static std::vector<double> GradientConjugate(std::vector<double>& x);
  
    
    
  }; 

}

