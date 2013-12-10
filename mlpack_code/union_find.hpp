
#ifndef UNION_FIND_HPP_
#define UNION_FIND_HPP_

#include <vector>
#include "data.hpp"

namespace bmst {

class UnionFind {

public:
  
  UnionFind(size_t num_points);
  
  ~UnionFind();
  
  size_t Find(const size_t& i);
  
  void Union(const size_t& i, const size_t& j);
  
  
private:

  std::vector<size_t> parents_;
  std::vector<size_t> ranks_;
  
}; // class

}


#endif
