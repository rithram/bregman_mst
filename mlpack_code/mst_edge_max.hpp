#ifndef MST_EDGE_MAX_HPP_
#define MST_EDGE_MAX_HPP_

namespace bmst {

  // A policy class which computes the length of an edge (x,y) as 
  // max(d_f(x,y), d_f(y,x)) for the given divergence f
  template<typename T, class TBregmanDiv>
  class MstMaxEdge {
  
    typedef BregmanBall<T, TBregmanDiv> BoundType;
  
  public:
    
    static double EdgeWeight(const Point<T>& x, const Point<T>&y);
  
    static bool CanPrune(const BoundType& query_bound, const BoundType& ref_bound);
                           
    static bool CanPrune(const Point<T>& query, const BoundType& ref_bound, double candidate_dist);                         
  
  }; // class



} // namespace




#endif