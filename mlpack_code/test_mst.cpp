
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
  std::uniform_real_distribution<double> randu(0, 10);
  
  std::vector<std::vector<double> > data_points;

  int num_data = 1000;
  int num_features = 5;

  for (size_t i = 0; i < num_data; i++) {
    std::vector<double> point;
    for (size_t j = 0; j < num_features; j++) {
      point.push_back(randu(generator));
    }
    data_points.push_back(point);
  }

  Table<double> data(data_points);
  
  
  //////////////////////////////////
  
  typedef KLDivergence<double> DivType;
  typedef BregmanBallTree<double, DivType, KMeansSplitter<double, DivType> > TreeType;
  
  MinimumSpanningTree<double, MstMaxEdge<double, DivType>, TreeType> mst(data);

  mst.ComputeNaive(true);

  std::vector<Edge> naive_edges = mst.EdgeList();

  return 0;
  
} // main

