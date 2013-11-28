/**
 * @file bmst/mlpack_code/test_kmeans_splitter.cpp
 * @author Pari Ram (p.ram@gatech.edu)
 *
 * This file tests the functionalities available to the kmeans
 * partitioning function.
 */

#include <assert.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "data.hpp"
#include "KLDivergence.hpp"
#include "L2Divergence.hpp"
#include "kmeans_splitter.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  // test where you are given 10 points with only 3 distinct points 
  // and the kmeans initialization should be able to catch this if 
  // we do 3-means clustering
  std::cout << "Testing the k-means initialization with L2Div ..." << std::endl;
  {
    std::vector<double> p1 {3, 0, 1};
    std::vector<double> p2 {9, 1, 2};
    std::vector<double> p3 {5, 1, 9};

    std::vector<std::vector<double> > test_points;
    test_points.push_back(p1);
    test_points.push_back(p2);
    test_points.push_back(p3);
    test_points.push_back(p1);
    test_points.push_back(p3);
    test_points.push_back(p2);
    test_points.push_back(p3);
    test_points.push_back(p2);
    test_points.push_back(p1);
    test_points.push_back(p2);
    // true membership (0, 1, 2, 0, 2, 1, 2, 1, 0, 1)

    const bmst::Table<double> test_table(test_points);
    bmst::KMeansSplitter<double, bmst::L2Divergence<double> > kmeans_test(3, 0);
    std::vector<size_t> membership;
    std::vector<bmst::Point<double> > centers;
    std::vector<double> radii;
    kmeans_test.PartitionData(
        test_table, 
        0, 
        test_table.n_points(), 
        membership, 
        centers, 
        radii);

    assert(centers.size() == radii.size());
    assert(centers.size() == 3);
    for (size_t i = 0; i < centers.size(); i++) 
    {
      std::cout << "Center " << i + 1 << " [ ";
      for (size_t j = 0; j < centers[i].n_dims(); j++)
        std::cout << centers[i][j] << " ";
      std::cout << "], Radius: " << radii[i] << std::endl;
    }

    // all radii should be 0
    for (size_t j = 0; j < radii.size(); j++)
      assert(radii[j] == 0.0);

    // check that the memberships are as expected
    // if 2 points have 0 divergence, they should have same membership
    for (size_t i = 0; i < test_table.n_points(); i++)
      for (size_t j = 0; j < test_table.n_points(); j++)
        if (bmst::L2Divergence<double>::Divergence(test_table[i], test_table[j]) == 0)
          assert(membership[i] == membership[j]);

    // make sure that the centers are either p1, p2 or p3
    std::vector<bool> points_seen {false, false, false};
    for (size_t i = 0; i < centers.size(); i++) 
    {
      bool point_matched = true;
      assert(centers[i].n_dims() == p1.size());
      for (size_t j = 0; j < p1.size(); j++)
        if (centers[i][j] != p1[j])
          point_matched = false;

      if (point_matched)
      {
        points_seen[0] = true;
        continue;
      }

      point_matched = true;
      assert(centers[i].n_dims() == p2.size());
      for (size_t j = 0; j < p2.size(); j++)
        if (centers[i][j] != p2[j])
          point_matched = false;

      if (point_matched)
      {
        points_seen[1] = true;
        continue;
      }

      point_matched = true;
      assert(centers[i].n_dims() == p3.size());
      for (size_t j = 0; j < p3.size(); j++)
        if (centers[i][j] != p3[j])
          point_matched = false;

      if (point_matched)
        points_seen[2] = true;
    }

    for (size_t i = 0; i < points_seen.size(); i++)
      assert(points_seen[i]);

  }
  std::cout << "Testing the k-means initialization with L2Div ... DONE" << std::endl;


  std::cout << "Testing the k-means initialization with KLDiv ..." << std::endl;
  {
    std::vector<double> p1 {3, 2, 1};
    std::vector<double> p2 {9, 1, 2};
    std::vector<double> p3 {5, 1, 9};

    std::vector<std::vector<double> > test_points;
    test_points.push_back(p1);
    test_points.push_back(p2);
    test_points.push_back(p3);
    test_points.push_back(p1);
    test_points.push_back(p3);
    test_points.push_back(p2);
    test_points.push_back(p3);
    test_points.push_back(p2);
    test_points.push_back(p1);
    test_points.push_back(p2);
    // true membership (0, 1, 2, 0, 2, 1, 2, 1, 0, 1)

    const bmst::Table<double> test_table(test_points);
    bmst::KMeansSplitter<double, bmst::KLDivergence<double> > kmeans_test(3, 0);
    std::vector<size_t> membership;
    std::vector<bmst::Point<double> > centers;
    std::vector<double> radii;
    kmeans_test.PartitionData(
        test_table, 
        0, 
        test_table.n_points(), 
        membership, 
        centers, 
        radii);

    assert(centers.size() == radii.size());
    assert(centers.size() == 3);
    for (size_t i = 0; i < centers.size(); i++) 
    {
      std::cout << "Center " << i + 1 << " [ ";
      for (size_t j = 0; j < centers[i].n_dims(); j++)
        std::cout << centers[i][j] << " ";
      std::cout << "], Radius: " << radii[i] << std::endl;
    }

    // all radii should be 0
    for (size_t j = 0; j < radii.size(); j++)
      assert(radii[j] == 0.0);

    // check that the memberships are as expected
    // if 2 points have 0 divergence, they should have same membership
    for (size_t i = 0; i < test_table.n_points(); i++)
      for (size_t j = 0; j < test_table.n_points(); j++)
        if (bmst::KLDivergence<double>::Divergence(test_table[i], test_table[j]) == 0)
          assert(membership[i] == membership[j]);

    // make sure that the centers are either p1, p2 or p3
    std::vector<bool> points_seen {false, false, false};
    for (size_t i = 0; i < centers.size(); i++) 
    {
      bool point_matched = true;
      assert(centers[i].n_dims() == p1.size());
      for (size_t j = 0; j < p1.size(); j++)
        if (centers[i][j] != p1[j])
          point_matched = false;

      if (point_matched)
      {
        points_seen[0] = true;
        continue;
      }

      point_matched = true;
      assert(centers[i].n_dims() == p2.size());
      for (size_t j = 0; j < p2.size(); j++)
        if (centers[i][j] != p2[j])
          point_matched = false;

      if (point_matched)
      {
        points_seen[1] = true;
        continue;
      }

      point_matched = true;
      assert(centers[i].n_dims() == p3.size());
      for (size_t j = 0; j < p3.size(); j++)
        if (centers[i][j] != p3[j])
          point_matched = false;

      if (point_matched)
        points_seen[2] = true;
    }

    for (size_t i = 0; i < points_seen.size(); i++)
      assert(points_seen[i]);

  }
  std::cout << "Testing the k-means initialization with KLDiv ... DONE" <<
    std::endl;

  // Test where there are 15 points with 6 distinct points and the 
  // 3 means clustering is done (where there are 3 well separated 
  // clusters) and the clustering membership should find these three 
  // clusters
  std::cout << "Testing the k-means clustering with L2Div ..." << std::endl;
  {
    bmst::Point<double> p11(std::vector<double>({3, 3}));
    bmst::Point<double> p12(std::vector<double>({3.1, 2.9}));
    bmst::Point<double> p21(std::vector<double>({9, 1}));
    bmst::Point<double> p22(std::vector<double>({8.8, 1.2}));
    bmst::Point<double> p31(std::vector<double>({1, 10}));
    bmst::Point<double> p32(std::vector<double>({1.1, 9}));

    std::vector<bmst::Point<double> > test_points;
    test_points.push_back(p11);// 0
    test_points.push_back(p21);// 1
    test_points.push_back(p31);// 2
    test_points.push_back(p12);// 0
    test_points.push_back(p31);// 2
    test_points.push_back(p22);// 1
    test_points.push_back(p31);// 2
    test_points.push_back(p21);// 1
    test_points.push_back(p12);// 0
    test_points.push_back(p22);// 1
    test_points.push_back(p21);// 1
    test_points.push_back(p32);// 2
    test_points.push_back(p22);// 1
    test_points.push_back(p11);// 0
    test_points.push_back(p12);// 0

    std::vector<bmst::Point<double> > true_centers;
    bmst::Point<double> center0 = (2 * p11 + 3 * p12) / 5;
    std::cout << "True center 1 [ " << center0[0] << " " << center0[1] <<
      " ]" << std::endl;
    bmst::Point<double> center1 = (3 * p21 + 3 * p22) / 6;
    std::cout << "True center 2 [ " << center1[0] << " " << center1[1] <<
      " ]" << std::endl;
    bmst::Point<double> center2 = (3 * p31 + p32) / 4;
    std::cout << "True center 3 [ " << center2[0] << " " << center2[1] <<
      " ]" << std::endl;

    const bmst::Table<double> test_table(test_points);
    bmst::KMeansSplitter<double, bmst::L2Divergence<double> > kmeans_test(3);
    std::vector<size_t> membership;
    std::vector<bmst::Point<double> > centers;
    std::vector<double> radii;
    kmeans_test.PartitionData(
        test_table, 
        0, 
        test_table.n_points(), 
        membership, 
        centers, 
        radii);

    assert(centers.size() == radii.size());
    assert(centers.size() == 3);
    for (size_t i = 0; i < centers.size(); i++) 
    {
      std::cout << "Center " << i + 1 << " [ ";
      for (size_t j = 0; j < centers[i].n_dims(); j++)
        std::cout << centers[i][j] << " ";
      std::cout << "], Radius: " << radii[i] << std::endl;
    }

    // make sure that the centers are as expected
    std::vector<bool> centers_seen {false, false, false};
    for (size_t i = 0; i < centers.size(); i++) 
    {
      bool center_matched = true;
      assert(centers[i].n_dims() == center0.n_dims());
      for (size_t j = 0; j < center0.n_dims(); j++)
      {
        if (fabs(centers[i][j] - center0[j]) > 0 and 
            fabs(centers[i][j] - center0[j]) < 1e-10)
          std::cout << "Center[" << i + 1 << "][" << j + 1 << "] diff " << 
            fabs(centers[i][j] - center0[j]) << std::endl;

        if (fabs(centers[i][j] - center0[j]) > 1e-10)
          center_matched = false;
      }

      if (center_matched)
      {
        centers_seen[0] = true;
        continue;
      }

      center_matched = true;
      assert(centers[i].n_dims() == center1.n_dims());
      for (size_t j = 0; j < center1.n_dims(); j++)
      {
        if (fabs(centers[i][j] - center1[j]) > 0 and 
            fabs(centers[i][j] - center0[j]) < 1e-10)
          std::cout << "Center[" << i + 1 << "][" << j + 1 << "] diff " << 
            fabs(centers[i][j] - center1[j]) << std::endl;

        if (fabs(centers[i][j] - center1[j]) > 1e-10)
          center_matched = false;
      }

      if (center_matched)
      {
        centers_seen[1] = true;
        continue;
      }

      center_matched = true;
      assert(centers[i].n_dims() == center2.n_dims());
      for (size_t j = 0; j < center2.n_dims(); j++)
      {
        if (fabs(centers[i][j] - center2[j]) > 0 and 
            fabs(centers[i][j] - center0[j]) < 1e-10)
          std::cout << "Center[" << i + 1 << "][" << j + 1 << "] diff " << 
            fabs(centers[i][j] - center2[j]) << std::endl;

        if (fabs(centers[i][j] - center2[j]) > 1e-10)
          center_matched = false;
      }

      if (center_matched)
        centers_seen[2] = true;

      if (not center_matched)
        std::cout << "Center " << i + 1 << " not matched !" << std::endl;
    }

    for (size_t i = 0; i < centers_seen.size(); i++)
      assert(centers_seen[i]);
  }
  std::cout << "Testing the k-means clustering with L2Div ... DONE" << 
    std::endl;

  // Test if the kmeans clustering works for a large set of random points
  // We will randomly select 1000 points in 5 dimensions and do a 
  // 2-means clustering
  std::cout << "Testing kmeans clustering on 1000 points with KLDiv ... " << 
    std::endl;
  {
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_real_distribution<double> urand(1e-10, 1);
    std::vector<bmst::Point<double> > point_set;
    for (size_t i = 0; i < 1000; i++)
    {
      std::vector<double> rand_vec;
      for (size_t j = 0; j < 5; j++)
        rand_vec.push_back(urand(gen));

      point_set.push_back(bmst::Point<double>(rand_vec));
    }
    const bmst::Table<double> rand_table(point_set);
    std::cout << "Clustering " << rand_table.n_points() << " points in " <<
      rand_table[0].n_dims() << " dimensions each .. " << std::endl;

    bmst::KMeansSplitter<double, bmst::KLDivergence<double> > kmeans_test(2);
    std::vector<size_t> membership;
    std::vector<bmst::Point<double> > centers;
    std::vector<double> radii;
    kmeans_test.PartitionData(
        rand_table, 
        0, 
        rand_table.n_points(), 
        membership, 
        centers, 
        radii);

    assert(membership.size() == rand_table.n_points());
    assert(centers.size() == radii.size());
    assert(centers.size() == 2);
    for (size_t i = 0; i < centers.size(); i++) 
    {
      std::cout << "Center " << i + 1 << " [ ";
      for (size_t j = 0; j < centers[i].n_dims(); j++)
        std::cout << centers[i][j] << " ";
      std::cout << "], Radius: " << radii[i] << std::endl;
    }

    std::vector<size_t> cluster_counts(2, 0);
    for (size_t i = 0; i < membership.size(); i++)
      cluster_counts[membership[i]]++;

    for (size_t i = 0; i < cluster_counts.size(); i++) 
      std::cout << "Cluster " << i + 1 << " has " << cluster_counts[i] << 
        " points." << std::endl;
  }
  std::cout << "Testing kmeans clustering on 1000 points with KLDiv ... " << 
    "DONE" << std::endl;


  return 0;
}
