
#include "left_nn_search.hpp"


using namespace bmst;

int main(int argc, char* argv[]) 
{
  
  Table<double> references("references.txt");
  Table<double> queries("queries.txt");

  // build the search class
  int leaf_size = 1;
  double min_ball_width = 1e-10;
  
  LeftNNSearch<Table<double>, Point<double>, KLDivergence<double> > searcher(references, leaf_size, min_ball_width);

  std::vector<size_t> neighbors(queries.n_points());

  for (int q = 0; q < queries.n_points(); q++)
  {
    
    neighbors[q] = searcher.ComputeNeighbor(queries[q]);
    
  } // loop over queries
  
  
  return 0;
  
}
