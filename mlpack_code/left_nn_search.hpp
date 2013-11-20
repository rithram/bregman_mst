

namespace bbtree {

  template<typename TDataType, class TBregmanDiv>
  class LeftNNSearch<TDataType, TBregmanDiv> {
  
  public:
    
    LeftNNSearch(std::vector<TDataType>& data, int leaf_size, double min_ball_width);
    
    ~LeftNNSearch();
    
    void Compute();
    
    
  private:
  
    std::vector<TDataType> data_;
    
    BregmanBallTree<TDataType, TBregmanDiv> tree_;

    int leaf_size_;
    double min_ball_width_;
    
    int neighbor_index_;
    double neighbor_distance_;
  
    std::vector<int> old_from_new_indices_;
    
    // functions
    
    void SearchNode_(BregmanBallTree<TDataType, TBregmanDiv>* node);
  
  }; // class

} // namespace