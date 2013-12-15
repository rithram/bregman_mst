
#include "minimum_spanning_tree.hpp"
#include "mst_edge_max.hpp"
#include "KLDivergence.hpp"
#include "L2Divergence.hpp"
#include "bregman_ball_tree.hpp"
#include "kmeans_splitter.hpp"

using namespace bmst;


int main(int argc, char* argv[])
{

  std::default_random_engine generator(time(NULL));
  std::uniform_real_distribution<double> randu(1, 10);
  
  std::vector<std::vector<double> > data_points;

  int num_data = 10;
  int num_features = 1;

  for (size_t i = 0; i < num_data; i++) {
    std::vector<double> point;
    for (size_t j = 0; j < num_features; j++) {
      point.push_back(randu(generator));
    }
    data_points.push_back(point);
  }

  Table<double> data(data_points);
  
  
  //////////////////////////////////
  
  typedef L2Divergence<double> DivType;
  typedef BregmanBallTree<double, DivType, KMeansSplitter<double, DivType> > TreeType;
  
  MinimumSpanningTree<double, MstMaxEdge<double, DivType>, TreeType> naive_false_mst(data, 1000);
  MinimumSpanningTree<double, MstMaxEdge<double, DivType>, TreeType> naive_true_mst(data, 1000);

  std::cout << "Comparing naive MST constructions\n";

  naive_false_mst.ComputeNaive(false);
  naive_true_mst.ComputeNaive(true);

  std::vector<Edge> naive_false_edges = naive_false_mst.EdgeList();
  std::vector<Edge> naive_true_edges = naive_true_mst.EdgeList();

  assert(naive_false_edges.size() == naive_true_edges.size());

  for (size_t i = 0; i < naive_false_edges.size(); i++)
  {
    
    assert((naive_false_edges[i].u == naive_true_edges[i].u and naive_false_edges[i].v == naive_true_edges[i].v) 
      or (naive_false_edges[i].u == naive_true_edges[i].v and naive_false_edges[i].v == naive_true_edges[i].u));
    
    assert(fabs(naive_false_edges[i].weight - naive_true_edges[i].weight) < 1e-6);
    
  }
  
  std::cout << "Naive constructions pass.\n\n";
  
  std::cout << "Testing Single-Tree Boruvka algorithm\n";
  
  MinimumSpanningTree<double, MstMaxEdge<double, DivType>, TreeType> single_mst(data);
  
  single_mst.ComputeSTB();

  std::vector<Edge> single_edges = single_mst.EdgeList();
  
  assert(single_edges.size() == naive_true_edges.size());

  for (size_t i = 0; i < naive_false_edges.size(); i++)
  {
    
    assert((single_edges[i].u == naive_true_edges[i].u and single_edges[i].v == naive_true_edges[i].v) 
      or (single_edges[i].u == naive_true_edges[i].v and single_edges[i].v == naive_true_edges[i].u));
    
    assert(fabs(single_edges[i].weight - single_edges[i].weight) < 1e-6);
    
  }
  
  std::cout << "Single-Tree Boruvka passes.\n";
  
  return 0;
  
} // main

