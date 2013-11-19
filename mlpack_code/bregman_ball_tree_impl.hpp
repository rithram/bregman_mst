

namespace bbtree {

  
  BregmanBallTree<TDataType, TBregmanDiv>BregmanBallTree(std::vector<TDataType>& data,
                                                         int begin, int count)
  :
  begin_(begin),
  end_(begin+count),
  count_(count),
  center_(),
  radius_(),
  left_(NULL),
  right_(NULL),
  bounding_ball_()
  {
    
    // update this to use right syntax for std::vector
    bounding_ball_ |= data(begin:end);
    
    center_ = bounding_ball_.center();
    radius_ = bounding_ball_.radius();
            
  }
  
  BregmanBallTree<TDataType, TBregmanDiv>~BregmanBallTree()
  {
    if (!IsLeaf())
    {
      ~left_();
      ~right_();
    }
    
    // free the bounding ball also, if necessary
    
  }
  

  bool BregmanBallTree<TDataType, TBregmanDiv>IsLeaf() const
  {
    if (left_) 
    {
      return true;
    }
    else {
      return false;
    }
  }
  
  BregmanBallTree<TDataStorage, TDataType>* BregmanBallTree<TDataType, TBregmanDiv>Left() const
  {
    return left_;
  }
  
  BregmanBallTree<TDataStorage, TDataType>* BregmanBallTree<TDataType, TBregmanDiv>Right() const
  {
    return right_;
  }
  
  int BregmanBallTree<TDataType, TBregmanDiv>Begin() const
  {
    return begin_;
  }

  int BregmanBallTree<TDataType, TBregmanDiv>End() const
  {
    return end_;
  }

  int BregmanBallTree<TDataType, TBregmanDiv>Count() const
  {
    return end_ - begin_;
  }
  
  int BregmanBallTree<TDataType, TBregmanDiv>Center() const
  {
    return center_;
  }
  
  double BregmanBallTree<TDataType, TBregmanDiv>Radius() const
  {
    return radius_;
  }
  
  BregmanBall<TBregmanDiv>& BregmanBallTree<TDataType, TBregmanDiv>Bound() const
  {
    return bounding_ball_;
  }
  
  
  
  template<typename TDataType, class TBregmanDiv>
  BregmanGallTree<TDataType, TBregmanDiv>* ConstructBBTree<TDataType, TBregmanDiv>(std::vector<TDataType>& data, int leaf_size, 
  double min_width)
  {
    
    BregmanBallTree<TDataType, TBregmanDiv>* tree = new BregmanBallTree<TDataType, TBregmanDiv>(0, data.size(), );
    
    SplitNode(tree, data, leaf_size, min_width);
    
    return tree;
    
  }

  template<typename TDataType, class TBregmanDiv>
  void SplitNode<TDataType, TBregmanDiv>(BregmanBallTree<TDataType, TBregmanDiv>* node, std::vector<TDataType>& data, 
                                         int leaf_size, double min_width)
  {

    // can we split at all?    
    if (node->count() <= leaf_size || node->radius() <= min_width) {
      return;
    }
    
    // Do 2-means clustering on the data subset represented by the node
    // This returns the number of points that belong to the left child
    int num_left = PartitionData(data(node.begin():node.end()));
    
    left = new BregmanBallTree<TDataType, TBregmanDiv>(data, node->begin(), num_left);
    right = new BregmanBallTree<TDataType, TBregmanDiv>(data, node->begin() + num_left, node->count() - num_left);
    
    SplitNode(left, data, leaf_size, min_width);
    SplitNode(right, data, leaf_size, min_width);
    
  } // SplitNode()
    
} // namespace