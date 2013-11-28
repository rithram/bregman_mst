/**
 *
 */

#ifndef BMST_KMEANS_SPLITTER_IMPL_HPP_
#define BMST_KMEANS_SPLITTER_IMPL_HPP_

#include <time.h>

#include <random>
#include <set>

#include "kmeans_splitter.hpp"

namespace bmst {

template<typename T, class TBregmanDiv>
KMeansSplitter<T, TBregmanDiv>::KMeansSplitter(
    const size_t k, const size_t max_iters) :
  k_(k),
  max_iterations_(max_iters)
{}

template<typename T, class TBregmanDiv>
void KMeansSplitter<T, TBregmanDiv>::PartitionData(
    const Table<T>& data,
    const size_t begin_index,
    const size_t end_index,
    std::vector<size_t>& membership,
    std::vector<Point<T> >& centers,
    std::vector<double>& radii) 
{
  // initialize the centers and such
  membership.resize(end_index - begin_index);
  size_t n_dims = data[begin_index].n_dims();
  // pick k_ random points and make them centers
  centers.resize(0);
  // first a random point to be the first center
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_int_distribution<size_t> urand(begin_index, end_index - 1);
  std::set<size_t> points_already_picked;
  size_t first_point = urand(gen);
  centers.push_back(data[first_point]);
  points_already_picked.insert(first_point);
  std::vector<double> div_to_closest_mean(
      end_index - begin_index, std::numeric_limits<double>::max());
  div_to_closest_mean[first_point] = 0;
  for (size_t j = 1; j < k_; j++) 
  {
    // compute distances to the previous center
    // and also keep track of the point with the largest div
    double max_div_to_closest_center = 0;
    size_t max_index = end_index;
    for (size_t i = begin_index; i < end_index; i++) 
    {
      double div_to_center = TBregmanDiv::Divergence(data[i], centers[j - 1]);
      if (div_to_center < div_to_closest_mean[i - begin_index])
      {
        div_to_closest_mean[i - begin_index] = div_to_center;
      }
      if (div_to_closest_mean[i - begin_index] > max_div_to_closest_center)
      {
        max_div_to_closest_center = div_to_closest_mean[i - begin_index];
        max_index = i;
      }
    }
    assert(max_index >= begin_index);

    if (max_index == end_index)
    {
      assert(max_div_to_closest_center == 0);
      std::cout << "[ERROR] Cannot obtain " << k_ << " distinct points " << 
        "to initialize " << k_ << " clusters. Can only get " << j << 
        " points." << std::endl;
      exit(1);
    }

    if (points_already_picked.find(max_index) == points_already_picked.end())
    {
      centers.push_back(data[max_index]);
      points_already_picked.insert(max_index);
      div_to_closest_mean[max_index - begin_index] = 0;
    }
    else
    {
      std::cout << "[ERROR] This should never happen because we should "
        "never pick the same point twice in this initialization scheme." <<
        std::endl;
      exit(1);
    }
  }

  // the re-assignment loop
  std::vector<size_t> old_membership;
  bool converged = false;
  std::vector<double> cluster_counts(k_);
  size_t num_iters = 0;
  double kmeans_obj;
  do
  {
    // assign the points to the centers
    if (membership.size() == 0)
      membership.resize(end_index - begin_index);

    kmeans_obj = 0;
    for (size_t i = begin_index; i < end_index; i++)
    {
      double min_div = std::numeric_limits<double>::max();
      size_t min_index = k_;
      for (size_t j = 0; j < k_; j++)
      {
        double div_to_center = TBregmanDiv::Divergence(data[i], centers[j]);
        if (div_to_center < min_div) 
        {
          min_div = div_to_center;
          min_index = j;
        }
      }
      assert(min_index < k_);
      membership[i - begin_index] = min_index;
      kmeans_obj += min_div;
    }
    // std::cout << "Obj: " << kmeans_obj << " @ iter " << 
    //   num_iters << std::endl;

    if (kmeans_obj == 0)
      converged = true;

    // compare new to old
    if (old_membership.size() > 0) 
    {
      bool all_equal = true;
      for (size_t i = 0; i < membership.size(); i++)
      {
        if (old_membership[i] != membership[i])
        {
          all_equal = false;
          break;
        }
      }
      // all new same all old
      if (all_equal) 
      {
        converged = true;
        break;
      }
    }

    // compute the new means for the assignment
    for (size_t j = 0; j < k_; j++)
      centers[j].zeros();
    cluster_counts.assign(k_, 0);

    for (size_t i = begin_index; i < end_index; i++)
    {
      centers[membership[i - begin_index]] += data[i];
      cluster_counts[membership[i - begin_index]]++;
    }

    for (size_t j = 0; j < k_; j++)
      if (cluster_counts[j] > 0)
        centers[j] /= cluster_counts[j];

    old_membership.swap(membership);     

  } while (num_iters++ < max_iterations_);

  if (converged and old_membership.size() > 0 and membership.size() > 0)
    for (size_t i = 0; i < membership.size(); i++)
      assert(membership[i] == old_membership[i]);
  else
  {
    if (not converged)
    {
      std::cout << "[WARNING] The k-means algorithm did not converge "
        "with " << max_iterations_ << " iterations." << std::endl;
      std::cout << "[WARNING] Using the most recent cluster assignment." <<
        std::endl;
    }
    if (membership.size() == 0)
      // just take the latest assignment
      membership.swap(old_membership);
  }
  // compute the radii for each of the centers
  radii.assign(k_, 0);
  for (size_t i = begin_index; i < end_index; i++)
  {
    size_t j = membership[i - begin_index];
    double div_to_center = TBregmanDiv::Divergence(data[i], centers[j]);
    if (div_to_center > radii[j])
      radii[j] = div_to_center;
  }
     
  return;
} // PartitionData

} // namespace

#endif
