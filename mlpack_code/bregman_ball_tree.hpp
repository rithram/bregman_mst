//
//  bregman_ball_tree.h
//  
//
//  Created by march on 11/12/13.
//
//

#ifndef ____bregman_ball_tree__
#define ____bregman_ball_tree__

#include "bregman_ball.hpp"

namespace bbtree {
  
  template <class TDataStorage, typename TDataType, class TBregmanDiv>
  class BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>
  {
    
  private:
    
    // Indices into the data set
    int begin_;
    int end_;
    int count_;
    
    // index of the center of the ball
    int center_;
    
    double radius_;
    
    BregmanBallTree<TDataStorage, TDataType>* left_;
    BregmanBallTree<TDataStorage, TDataType>* right_;
    
    
    BregmanBall<TBregmanDiv> bounding_ball_;
    
    
  public:
    
    BregmanBallTree();
    
    ~BregmanBallTree();
    

    bool IsLeaf() const;
    
    BregmanBallTree<TDataStorage, TDataType>* Left() const;
    BregmanBallTree<TDataStorage, TDataType>* Right() const;
    
    int Begin() const;
    int End() const;
    int Count() const;
    
    int Center() const;
    TDataType Center() const;
    
    double Radius() const;
    
    BregmanBall<TBregmanDiv>& Bound() const;
    
    
      
      
    
    
    
    
  }; // class
  
} // namespace

#include "bregman_ball_tree_impl.hpp"







#endif /* defined(____bregman_ball_tree__) */
