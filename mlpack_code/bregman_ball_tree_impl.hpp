

namespace bbtree {

  
  BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>BregmanBallTree();
  
  BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>~BregmanBallTree();
  

  bool BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>IsLeaf() const;
  
  BregmanBallTree<TDataStorage, TDataType>* BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>Left() const;
  BregmanBallTree<TDataStorage, TDataType>* BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>Right() const;
  
  int BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>Begin() const;
  int BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>End() const;
  int BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>Count() const;
  
  int BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>Center() const;
  TDataType BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>Center() const;
  
  double BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>Radius() const
  {
    return radius_;
  }
  
  BregmanBall<TBregmanDiv>& BregmanBallTree<TDataStorage, TDataType, TBregmanDiv>Bound() const
  {
    return bounding_ball_;
  }
  


}