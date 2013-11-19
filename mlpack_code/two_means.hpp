/**
 * 2-means clustering for constructing Bregman ball trees.
 */
 
 namespace bbtree {


   template<typename TDataType, class TBregmanDiv>
   int PartitionData<TDataType, TBregmanDiv>(std::vector<TDataType>& data);
 
 
 }
 
#include "two_means_impl.hpp"
 