#ifndef MST_EDGE_MAX_IMPL_HPP_
#define MST_EDGE_MAX_IMPL_HPP_

namespace bmst {

  template<typename T, class TBregmanDiv>
  static double MstEdgeMax<T, TBregmanDiv>::EdgeWeight(const Point<T>& x, const Point<T>&y)
  {
  
    double xy = TBregmanDiv::Divergence(x,y);
    double yx = TBregmanDiv::Divergence(y,x);
    
    return std::max(xy, yx);
    
  }

  template<typename T, class TBregmanDiv>
  static bool MstEdgeMax<T, TBregmanDiv>::CanPrune(const BoundType& query_bound,
                                                   const BoundType& ref_bound)
  { 

    return false;

  }

  template<typename T, class TBregmanDiv>
  static bool MstEdgeMax<T, TBregmanDiv>::CanPrune(const Point<T>& query,
                                                    const BoundType& ref_bound, 
                                                    double candidate_dist)
  {
    return false;
  }



}


#endif