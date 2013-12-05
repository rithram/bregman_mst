
#include "union_find.hpp"

using namespace bmst;

UnionFind::UnionFind(size_t num_points)
  :
parents_(num_points),
ranks_(num_points)
{
  
  for (size_t i = 0; i < num_points; i++)
  {
    parents_[i] = i;
    ranks_[i] = 0;
  }
}

UnionFind::~UnionFind()
{}

size_t UnionFind::Find(const size_t& i)
{
  if (parents_[i] == i) 
  {
    return i;
  }
  else 
  {
    parents_[i] = Find(parents_[i]);
    return parents_[i];
  }
}

void UnionFind::Union(const size_t& i, const size_t& j)
{
  
  size_t root_i = Find(i);
  size_t root_j = Find(j);
  
  if (root_i == root_j) {
    return; // they're already linked
  }
  else if (ranks_[root_i] == ranks_[root_j]) 
  {
    parents_[root_j] = root_i;
    ranks_[root_i]++;
  }
  else if (ranks_[root_i] < ranks_[root_j]) {
    parents_[root_i] = root_j;
  }
  else {
    parents_[root_j] = root_i;
  }
  
}


