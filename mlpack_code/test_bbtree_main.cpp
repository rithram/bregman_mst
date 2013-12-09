
/**
 * @file bmst/mlpack_code/test_bbtree_main.cpp
 * @author Pari Ram (p.ram@gatech.edu)
 *
 * This file tests the construction and behavior of the bbtree class
 */

#include <assert.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include "data.hpp"
#include "KLDivergence.hpp"
#include "L2Divergence.hpp"
#include "kmeans_splitter.hpp"
#include "bregman_ball_tree.hpp"

template <typename T, class TNode, class TBregmanDiv>
void TestTreeNode(const bmst::Table<T>& table, const TNode* node);

int main(int argc, char* argv[])
{
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<double> randu(1e-10, 10);
  std::cout << "================================================" << std::endl;
  std::cout << "Testing the bbtree with L2Div ... " << std::endl;
  {
    // make a 100 x 5 dataset
    std::vector<bmst::Point<double> > point_set;
    for (size_t i = 0; i < 1000; i++)
    {
      std::vector<double> rand_vec;
      for (size_t j = 0; j < 50; j++)
        rand_vec.push_back(randu(gen));

      point_set.push_back(bmst::Point<double>(rand_vec));
    }
    bmst::Table<double> rand_table(point_set);
    std::cout << "Indexing " << rand_table.n_points() << " points in " <<
      rand_table[0].n_dims() << " dimensions each .. " << std::endl;

    typedef bmst::L2Divergence<double> TBregmanDiv;
    typedef bmst::KMeansSplitter<double, TBregmanDiv> TSplitter;
    typedef bmst::BregmanBallTree<double, TBregmanDiv, TSplitter> BBTree;

    std::vector<size_t> old_from_new;
    BBTree* test_bbtree = new BBTree(rand_table, old_from_new, 5);

    // test the stats of each node
    std::queue<BBTree*> node_queue;
    node_queue.push(test_bbtree);
    while (not node_queue.empty())
    {
      BBTree* current_node = node_queue.front();
      TestTreeNode<double, BBTree, TBregmanDiv>(rand_table, current_node);
      node_queue.pop();
      if (not current_node->IsLeaf())
      {
        node_queue.push(current_node->Left());
        node_queue.push(current_node->Right());
      }
    }

    delete test_bbtree;
  }
  std::cout << "Testing the bbtree with L2Div ... DONE" << std::endl;

  std::cout << "================================================" << std::endl;
  std::cout << "Testing the bbtree with KLDiv ... " << std::endl;
  {
    // make a 1000 x 5 dataset
    std::vector<bmst::Point<double> > point_set;
    for (size_t i = 0; i < 1000; i++)
    {
      std::vector<double> rand_vec;
      for (size_t j = 0; j < 100; j++)
        rand_vec.push_back(randu(gen));

      point_set.push_back(bmst::Point<double>(rand_vec));
    }
    bmst::Table<double> rand_table(point_set);
    std::cout << "Indexing " << rand_table.n_points() << " points in " <<
      rand_table[0].n_dims() << " dimensions each .. " << std::endl;

    typedef bmst::KLDivergence<double> TBregmanDiv;
    typedef bmst::KMeansSplitter<double, TBregmanDiv> TSplitter;
    typedef bmst::BregmanBallTree<double, TBregmanDiv, TSplitter> BBTree;

    std::vector<size_t> old_from_new;
    BBTree* test_bbtree = new BBTree(rand_table, old_from_new, 5);

    // test the stats of each node
    std::queue<BBTree*> node_queue;
    node_queue.push(test_bbtree);
    while (not node_queue.empty())
    {
      BBTree* current_node = node_queue.front();
      TestTreeNode<double, BBTree, TBregmanDiv>(rand_table, current_node);
      node_queue.pop();
      if (not current_node->IsLeaf())
      {
        node_queue.push(current_node->Left());
        node_queue.push(current_node->Right());
      }
    }

    delete test_bbtree;
  }
  std::cout << "Testing the bbtree with KLDiv ... DONE" << std::endl;
  std::cout << "================================================" << std::endl;

  std::cout << "[TESTS-TO-BE-ADDED] We need to add tests for 'CentroidPrimes' and "
    "for the left center and left radius" << std::endl;

  return 0;
}

template <typename T, class TNode, class TBregmanDiv>
void TestTreeNode(const bmst::Table<T>& table, const TNode* node)
{
  // check the node center
  bmst::Point<T> center;
  center.zeros(table[0].n_dims());

  for (size_t i = node->Begin(); i < node->End(); i++)
    center += table[i];

  center /= (T) node->Count();

  assert(center.n_dims() == node->RCenter().n_dims());
  const bmst::Point<T> center_diff = center - node->RCenter();
  double center_diff_sq_norm = bmst::Dot(center_diff, center_diff);
  if (center_diff_sq_norm > 1e-10)
  {
    std::cout << "[TEST-FAIL] ||true center - supposed center||^2 = " << 
      center_diff_sq_norm << std::endl;
    exit(1);
  }

  // check the node radius
  double radius = 0;
  for (size_t i = node->Begin(); i < node->End(); i++)
  {
    double div_to_center = TBregmanDiv::Divergence(table[i], center);
    if (div_to_center > radius)
      radius = div_to_center;
  }

  if (fabs(radius - node->RRadius()) > 1e-10)
  {
    std::cout << "[TEST-FAIL] |true radius - supposed radius| = " << 
      fabs(radius - node->RRadius()) << std::endl;
    exit(1);
  }

  return;
}
