/**
 * @file 
 * 2-means clustering for constructing Bregman ball trees.
 */

#ifndef BMST_KMEANS_SPLITTER_HPP_
#define BMST_KMEANS_SPLITTER_HPP_

#include <vector>

#include "data.hpp"

namespace bmst {

template <typename T, class TBregmanDiv>
class KMeansSplitter
{
private:
  size_t k_;
  size_t max_iterations_;

public:
  KMeansSplitter(const size_t k = 2, const size_t max_iters = 10000);

  void PartitionData(
      const Table<T>& data,
      const size_t begin_index,
      const size_t end_index,
      std::vector<size_t>& membership,
      std::vector<Point<T> >& centers,
      std::vector<double>& radii);
}; // class KMeansSplitter
 
} // namespace
#include "kmeans_splitter_impl.hpp"
 
#endif
