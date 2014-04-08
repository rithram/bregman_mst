/**
 * @file bmst/mlpack_code/util.hpp
 * @author Pari Ram (p.ram@gatech.edu)
 *
 * Some random utility functions
 */

#ifndef BMST_UTIL_HPP_
#define BMST_UTIL_HPP_

#include <memory>

#include "data.hpp"

namespace bmst {
namespace util {

template <class TPoint>
bool PointHasZero(const TPoint& p);

template <typename T>
void SplitSet(
    const Table<T>& table, 
    const double split_ratio, 
    std::unique_ptr<Table<T> >& qset, 
    std::unique_ptr<Table<T> >& rset);

}; // namespace
}; // namespace

#endif

#include "util_impl.hpp"
