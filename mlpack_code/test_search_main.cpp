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

#include "left_nn_search.hpp"

using namespace std;

template <typename T, class Divergence>
void DoSearchAndCompareToNaive(bmst::Table<T>& rset, bmst::Table<T>& qset);

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
    DoSearchAndCompareToNaive<float, bmst::KLDivergence<float> >(*rset, *qset);
  else
  {  
    assert(chosen_divergence == "L2");
    DoSearchAndCompareToNaive<float, bmst::L2Divergence<float> >(*rset, *qset);
  }

  if (results_file != "")
  {
    cout << "The neighbors (and the corresponding divergences) are output in "
      "'" << results_file << "'" << endl;
  }

  return 0;
} // main


template <typename T, class Divergence>
void DoSearchAndCompareToNaive(bmst::Table<T>& rset, bmst::Table<T>& qset)
{
  bmst::LeftNNSearch<T, Divergence> searcher(rset, 1);
  std::vector<size_t> neighbors(qset.n_points());
  std::vector<size_t> naive_neighbors(qset.n_points());

  cout << "Testing search correctness ... ";

  for (size_t i = 0; i < qset.n_points(); i++) 
  {
    neighbors[i] = searcher.ComputeNeighbor(qset[i]);
    naive_neighbors[i] = searcher.ComputeNeighborNaive(qset[i]);
    assert(neighbors[i] == naive_neighbors[i]);
  }
  cout << "DONE " << endl;

  return;
}
