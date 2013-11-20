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
#include "two_means.hpp"

namespace bbtree {
  
  template <typename TDataType, class TBregmanDiv>
  class BregmanBallTree<TDataType, TBregmanDiv>
  {
    
  private:
    
    // Indices into the data set
    int begin_;
    int end_;
    int count_;
    
    // index of the center of the ball
    int center_;
    
    double radius_;
    
    BregmanBallTree<TDataType, TBregmanDiv>* left_;
    BregmanBallTree<TDataType, TBregmanDiv>* right_;
    
    
    BregmanBall<TBregmanDiv> bounding_ball_;
    
    
  public:
    
    BregmanBallTree(std::vector<TDataType>& data);
    
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
  
  template<typename TDataType, class TBregmanDiv>
  BregmanBallTree<TDataType, TBregmanDiv>* ConstructBBTree<TDataType, TBregmanDiv>(std::vector<TDataType>& data,
                                                                                   int leaf_size, double min_ball_width);
  
} // namespace

#include "bregman_ball_tree_impl.hpp"







#endif /* defined(____bregman_ball_tree__) */
