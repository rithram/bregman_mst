/**
 * @file bmst/mlpack_code/test_search_main.cpp
 * @author Pari Ram (p.ram@gatech.edu)
 *
 */

#include <assert.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "data.hpp"
#include "util.hpp"

#include "L2Divergence.hpp"
#include "KLDivergence.hpp"
#include "bregman_ball.hpp"
#include "left_nn_search.hpp"

using namespace std;

template <typename T, class Divergence>
void DoSearchAndCompareToNaive(
    bmst::Table<T>& rset, bmst::Table<T>& qset, const size_t leaf_size);

int main(int argc, char* argv[])
{
  namespace bpo = boost::program_options;

  // input command line options
  bpo::options_description opt_desc(
    "Options for performing searches with Bregman divergences");
  opt_desc.add_options()
    ("help", "Produce help message")
    ("rfile", bpo::value<string>(), 
     "The file containing the set of points (required)")
    ("qfile", bpo::value<string>(), 
     "A file containing a separate set of queries (optional)")
    ("divergence", bpo::value<string>(), 
     "The divergence to be used for the search (optional). Options are: \n"
     " L2 (default)\n"
     " KL\n")
    ("results", bpo::value<string>(), "The file in which to write the results")
    ("k", bpo::value<string>(),
     "The number of neighbors required for each query "
     "(optional, 'k' defaults to 1)")
    ("leaf_size", bpo::value<string>(),
     "The maximum number of points in any leaf of the tree "
     "(optional, 'leaf_size' defaults to 10)")
    ("split_ratio", bpo::value<string>(), "The ratio with which the dataset "
     "is split into query and reference sets (optional, defaults to 0.1 "
     "if the query set is not provided)");

  // read command line arguments
  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(argc, argv, opt_desc), vm);
  bpo::notify(vm);

  if (vm.count("help")) 
  {
    cout << opt_desc << endl;
    exit(0);
  }

  if (vm.count("rfile") == 0)
  {
    cout << "[ERROR] The --rfile option is required for specifying the set "
      "of points on which to perform the search on"  << endl;
    exit(1);
  }

  // Currently supported divergences
  set<string> divergences;
  divergences.insert("L2");
  divergences.insert("KL");

  if (vm.count("qfile") and vm.count("split_ratio")) 
    cout << "[WARNING] The query file is already provided, ignoring the "
      "--split_ratio provided." << endl;

  string rfile = vm["rfile"].as<string>();
  string qfile = vm.count("qfile") ? vm["qfile"].as<string>() : "";
  string chosen_divergence = 
    vm.count("divergence") ? vm["divergence"].as<string>() : "L2";
  string results_file = vm.count("results") ? vm["results"].as<string>() : "";
  size_t k = vm.count("k") ? atoi(vm["k"].as<string>().c_str()) : 1;
  double query_ref_split_ratio = vm.count("split_ratio") ? 
    atof(vm["split_ratio"].as<string>().c_str()) : 0.1;
  size_t leaf_size = vm.count("leaf_size") ? 
    atoi(vm["leaf_size"].as<string>().c_str()) : 10;

  if (divergences.find(chosen_divergence) == divergences.end())
  {
    // an unsupported divergence is selected
    cout << "[ERROR] " << chosen_divergence << 
      "-divergence is currently not supported" << endl;
    exit(1);
  }

  cout << "Reading in '" << rfile << "'" << endl;
  bmst::Table<float> data(rfile);
  std::unique_ptr<bmst::Table<float> > qset;
  std::unique_ptr<bmst::Table<float> > rset;

  if (qfile != "")
  {
    cout << "Reading in '" << qfile << "'" << endl;
    qset.reset(new bmst::Table<float>(qfile));
    rset.reset(new bmst::Table<float>(data));
  }
  else
  {
    // Splitting the dataset into a query and reference sets for the given ratio
    if (query_ref_split_ratio < 1.0) 
    {
      cout << "[INFO] Splitting the dataset randomly in the ratio " << 
        query_ref_split_ratio << " : " << 1.0 - query_ref_split_ratio << endl;
      bmst::util::SplitSet(data, query_ref_split_ratio, qset, rset);
    }
    else
    {
      cout << "[ERROR] The --split_ratio should be in the range (0, 1)" << endl;
      exit(0);
    }
  }

  cout << "Finding " << k << " neighbor(s) for " << qset->n_points() << 
    " queries from a set of " << rset->n_points() << 
    " points with respect to the " << chosen_divergence << "-divergence" << 
    endl;

  if (chosen_divergence == "KL")
    DoSearchAndCompareToNaive<float, bmst::KLDivergence<float> >(
        *rset, *qset, leaf_size);
  else
  {  
    assert(chosen_divergence == "L2");
    DoSearchAndCompareToNaive<float, bmst::L2Divergence<float> >(
        *rset, *qset, leaf_size);
  }

  if (results_file != "")
  {
    cout << "The neighbors (and the corresponding divergences) are output in "
      "'" << results_file << "'" << endl;
  }

  return 0;
} // main

template <typename T, class TDivergence>
void DoSearchAndCompareToNaive(
    bmst::Table<T>& rset, bmst::Table<T>& qset, const size_t leaf_size)
{
  //qset.make_non_zero(0.01);
  //rset.make_non_zero(0.02);

  cout << "[INFO] Indexing the reference set with leaves of maximum size " << 
    leaf_size << " ..." << endl;  
  typedef bmst::BregmanBall<T, TDivergence> TBBall;
  bmst::LeftNNSearch<T, TDivergence, TBBall> searcher(rset, leaf_size);
  cout << "[INFO] Reference set indexed" << endl;

  std::vector<size_t> neighbors(qset.n_points());
  std::vector<size_t> naive_neighbors(qset.n_points());

  cout << "[INFO] Testing search correctness ... ";
  size_t errors = 0;
  size_t num_queries_with_zero = 0;
  size_t total_bdiv_counter = 0;
  size_t total_grad_counter = 0;
  size_t total_grad_con_counter = 0;

  for (size_t i = 0; i < qset.n_points(); i++) 
  {
    TDivergence::bdiv_counter = 0;
    TDivergence::grad_counter = 0;
    TDivergence::grad_con_counter = 0;
    naive_neighbors[i] = searcher.ComputeNeighborNaive(qset[i]);
    if (naive_neighbors[i] == -1) {
      assert(bmst::util::PointHasZero(qset[i]));
      ++num_queries_with_zero;
    } else {
      assert(naive_neighbors[i] < rset.n_points());
      assert(TDivergence::bdiv_counter == rset.n_points());
      assert(TDivergence::grad_counter == 0);
      assert(TDivergence::grad_con_counter == 0);
    }
    TDivergence::bdiv_counter = 0;
    TDivergence::grad_counter = 0;
    TDivergence::grad_con_counter = 0;
    neighbors[i] = searcher.ComputeNeighbor(qset[i]);
    if (neighbors[i] == -1) {
      assert(bmst::util::PointHasZero(qset[i]));
    } else {
      assert(neighbors[i] < rset.n_points());
      total_bdiv_counter += TDivergence::bdiv_counter;
      total_grad_counter += TDivergence::grad_counter;
      total_grad_con_counter += TDivergence::grad_con_counter;
    }

    if (neighbors[i] != naive_neighbors[i]) 
    {
      if (neighbors[i] == -1 or naive_neighbors[i] == -1) {
        ++errors;
      } else {
        T naive_div = TDivergence::BDivergence(qset[i], rset[naive_neighbors[i]]);
        T div = TDivergence::BDivergence(qset[i], rset[neighbors[i]]);
        if (naive_div < div) ++errors;
      }
      // qset[i].print();
    }
  }
  cout << "DONE " << endl;
  if (errors > 0) 
    cout << "[ERROR] ";
  else
    cout << "[INFO] ";
  cout << errors << "/" << qset.n_points() << " errors" << endl;
  cout << "[INFO] " << num_queries_with_zero << "/" << qset.n_points() << 
    " queries with zero" << endl;
  cout << "[INFO] Naive comp: " << "D " << rset.n_points() * 
    (qset.n_points() - num_queries_with_zero) << endl;
  cout << "[INFO] Tree comp:  " << "D " << total_bdiv_counter << " G " << 
    total_grad_counter << " C " << total_grad_con_counter << endl;
  return;
}
