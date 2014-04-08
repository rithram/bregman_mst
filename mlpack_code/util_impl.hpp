/**
 * @file bmst/mlpack_code/util_impl.hpp
 * @author Pari Ram (p.ram@gatech.edu)
 *
 * Implementation of some random utility functions
 */

#ifndef BMST_UTIL_IMPL_HPP_
#define BMST_UTIL_IMPL_HPP_

#include <assert.h>
#include <time.h>

#include <algorithm>
#include <random>
#include <vector>

#include "util.hpp"

namespace bmst {
namespace util {

template <class TPoint>
bool PointHasZero(const TPoint& p) {
  bool has_zero = false;
  for (size_t i = 0; i < p.n_dims(); ++i) {
    if (p[i] == 0) {
      has_zero = true;
      break;
    }
  }
  return has_zero;
}

template <typename T>
void SplitSet(
    const Table<T>& table, 
    const double split_ratio, 
    std::unique_ptr<Table<T> >& qset, 
    std::unique_ptr<Table<T> >& rset)
{
  size_t num_queries = std::ceil(split_ratio * table.n_points());
  std::vector<size_t> qlist(num_queries, 0);
  std::vector<size_t> full_list(table.n_points() - num_queries, 1);
  full_list.insert(full_list.end(), qlist.begin(), qlist.end());
  assert(full_list.size() == table.n_points());

  // shuffle the list
  std::random_device rd;
  std::shuffle(
      full_list.begin(), full_list.end(), std::default_random_engine(rd()));

  std::vector<Point<T> > point_set;

  // Doing queries and references one by one to save memory
  // pick the queries
  for (size_t i = 0; i < table.n_points(); i++) 
  {
    if (full_list[i] == 0) 
      point_set.push_back(table[i]);
  }
  qset.reset(new Table<T>(point_set));
  point_set.clear();

  // pick the references
  for (size_t i = 0; i < table.n_points(); i++) 
  {
    if (full_list[i] == 1)
      point_set.push_back(table[i]);
  }
  rset.reset(new Table<T>(point_set));
  point_set.clear();

  return;
}

}; // namespace
}; // namespace

#endif
