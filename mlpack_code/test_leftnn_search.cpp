
#include "left_nn_search.hpp"
#include "KLDivergence.hpp"
#include "L2Divergence.hpp"

using namespace bmst;

int main(int argc, char* argv[]) 
{
  
  std::default_random_engine generator(time(NULL));
  std::uniform_real_distribution<double> randu(0, 10);
  
  std::vector<std::vector<double> > reference_points;
  std::vector<std::vector<double> > query_points;

  int num_references = 1000;
  int num_queries = 50;
  int num_features = 10;

  for (size_t i = 0; i < num_references; i++) {
    std::vector<double> ref_point;
    for (size_t j = 0; j < num_features; j++) {
      ref_point.push_back(randu(generator));
    }
    reference_points.push_back(ref_point);
  }

  for (size_t i = 0; i < num_queries; i++) {
    std::vector<double> query_point;
    for (size_t j = 0; j < num_features; j++) {
      query_point.push_back(randu(generator));
    }
    query_points.push_back(query_point);
  }

  Table<double> references(reference_points);
  Table<double> queries(query_points);

  // build the search class
  size_t leaf_size = 2;
  
  LeftNNSearch<double, KLDivergence<double> > searcher(references, leaf_size);

  std::vector<size_t> neighbors(queries.n_points());
  std::vector<size_t> naive_neighbors(queries.n_points());

  std::cout << "Testing KL Divergence Search.\n";

  for (int q = 0; q < queries.n_points(); q++)
  {
    
    neighbors[q] = searcher.ComputeNeighbor(queries[q]);

    naive_neighbors[q] = searcher.ComputeNeighborNaive(queries[q]);
    
    assert(neighbors[q] == naive_neighbors[q]);
    
  } // loop over queries

  std::cout << "KL Divergence tests PASSED.\n";

  leaf_size = 5;

  LeftNNSearch<double, L2Divergence<double> > searcher_l2(references, leaf_size);
  
  neighbors.clear();
  neighbors.resize(queries.n_points());
  naive_neighbors.clear();
  naive_neighbors.resize(queries.n_points());
  
  std::cout << "Testing L2 Divergence Search.\n";
  
  
  for (int q = 0; q < queries.n_points(); q++)
  {
    
    naive_neighbors[q] = searcher_l2.ComputeNeighborNaive(queries[q]);
    neighbors[q] = searcher_l2.ComputeNeighbor(queries[q]);
    
    assert(neighbors[q] == naive_neighbors[q]);
    
  }
  
  std::cout << "L2 Divergence tests PASSED.\n";
  
  
  return 0;
  
}
