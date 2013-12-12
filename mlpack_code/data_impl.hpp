/**
 * @file bregman_mst/mlpack_code/data_impl.hpp
 * @author Pari Ram (p.ram@gatech.edu)
 *
 * Implementation of the functions defined in data.hpp
 */

#ifndef BMST_DATA_IMPL_HPP_
#define BMST_DATA_IMPL_HPP_

#include <fstream> 
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "data.hpp"

namespace bmst
{ 

template <typename T>
Point<T>::Point() :
  values_(std::vector<T>(0)),
  n_dims_(0) 
{}

template <typename T>
Point<T>::Point(const std::vector<T>& point) :
  values_(point),
  n_dims_(point.size()) 
{}

template<typename T>
Point<T>::Point(const Point<T>& point)
{
  *this = point;
}

template<typename T>
T& Point<T>::operator[](const size_t i)
{
  if (i >= n_dims_) 
  {
    std::cout << "[ERROR] Point index out of range" << std::endl;
    exit(1);
  }
  return values_[i];
}

template<typename T>
const T& Point<T>::operator[](const size_t i) const
{
  if (i >= n_dims_) 
  {
    std::cout << "[ERROR] Point index out of range" << std::endl;
    exit(1);
  }
  return values_[i];
}

template<typename T>
Point<T>& Point<T>::operator=(const Point<T>& point)
{
  n_dims_ = point.n_dims_;
  values_ = point.values_;
  // values_.resize(n_dims_);
  // for (size_t i = 0; i < n_dims_; i++) 
  //   values_[i] = point.values_[i];

  return *this;
}

template<typename T>
Point<T>& Point<T>::operator+=(const Point<T>& point)
{
  if (n_dims_ == point.n_dims_) 
  {
    for (size_t i = 0; i < n_dims_; i++)
      values_[i] += point.values_[i];
  }
  else
  {
    std::cout << "[ERROR] Dimensions mismatch in point addition" <<
      std::endl;
    exit(1);
  }
  return *this;
}

template<typename T>
Point<T>& Point<T>::operator-=(const Point<T>& point)
{
  if (n_dims_ == point.n_dims_) 
  {
    for (size_t i = 0; i < n_dims_; i++)
      values_[i] -= point.values_[i];
  }
  else
  {
    std::cout << "[ERROR] Dimensions mismatch in point subtraction" <<
      std::endl;
    exit(1);
  }
  return *this;
}

template<typename T>
Point<T>& Point<T>::operator*=(const double scalar)
{
  for (size_t i = 0; i < n_dims_; i++)
    values_[i] *= scalar;

  return *this;
}

template<typename T>
Point<T>& Point<T>::operator/=(const double scalar)
{
  for (size_t i = 0; i < n_dims_; i++)
    values_[i] /= scalar;

  return *this;
}

template<typename T>
void Point<T>::zeros()
{
  values_.assign(n_dims_, 0);
}

template<typename T>
void Point<T>::zeros(const size_t num_dims)
{
  values_.assign(num_dims, 0);
  n_dims_ = num_dims;
}

template<typename T>
void Point<T>::ones()
{
  values_.assign(n_dims_, 1);
}

template<typename T>
void Point<T>::ones(const size_t num_dims)
{
  values_.assign(num_dims, 1);
  n_dims_ = num_dims;
}

template<typename T>
void Point<T>::print() const
{
  
  std::cout << "(";
  for (int i = 0; i < n_dims_ - 1; i++) 
  {
    std::cout << values_[i] << ", ";
  }
  std::cout << values_[n_dims_ - 1] << ")\n";
  
}

// Binary +
template<typename T>
Point<T> operator+(const Point<T>& lhs, const Point<T>& rhs)
{
  return Point<T>(lhs) += rhs;
}

// Binary -
template<typename T>
Point<T> operator-(const Point<T>& lhs, const Point<T>& rhs)
{
  return Point<T>(lhs) -= rhs;
}

// Binary *
template<typename T>
Point<T> operator*(const double lhs, const Point<T>& rhs)
{
  return Point<T>(rhs) *= lhs;
}

template<typename T>
Point<T> operator*(const Point<T>& lhs, const double rhs)
{
  return Point<T>(lhs) *= rhs;
}

// Binary /
template<typename T>
Point<T> operator/(const Point<T>& lhs, const double rhs)
{
  return Point<T>(lhs) /= rhs;
}

// Binary dot product
template<typename T>
double Dot(const Point<T>& a, const Point<T>& b)
{
  if (a.n_dims() != b.n_dims())
  {
    // fail here
    std::cout << "[ERROR] Dimension mismatch in dot product computation" <<
      std::endl;
    exit(1);
  }
  double dot_product = 0;
  // FIXME: Improve efficiency of this loop
  for (size_t i = 0; i < a.n_dims(); i++)
    dot_product += (a[i] * b[i]);
  return dot_product;
}

template<typename T>
Table<T>::Table() :
  points_(std::vector<Point<T> >(0)),
  n_points_(0)
{}

template<typename T>
Table<T>::Table(const std::vector<std::vector<T> >& points)
{
  for (size_t i = 0; i < points.size(); i++)
    points_.push_back(Point<T>(points[i]));
  n_points_ = points.size();
}

template<typename T>
Table<T>::Table(const std::vector<Point<T> >& points) :
  points_(points),
  n_points_(points.size())
{}

template<typename T>
Table<T>::Table(const Table<T>& table)
{
  *this = table;
}

template<typename T>
Table<T>::Table(const std::string& file_name) 
{
  points_.resize(0);
  n_points_ = 0;
  std::ifstream ifs;
  ifs.open(file_name.c_str(), std::ifstream::in);
  std::string line;
  boost::char_separator<char> sep(" ,");
  boost::tokenizer<boost::char_separator<char> > tok(line, sep);
  std::vector<std::string> line_pieces;
  std::vector<T> dims;
  size_t n_dims = 0;
  while (ifs.good()) {
    std::getline(ifs, line);
    tok.assign(line.begin(), line.end());
    line_pieces.assign(tok.begin(), tok.end());
    if (line_pieces.size() > 0) 
    {
      if (n_dims == 0) 
        n_dims = line_pieces.size();
      else 
      {
        if (n_dims != line_pieces.size()) 
        {
          std::cout << "[ERROR] Dimensionality of the points do not match." << 
            std::endl;
          exit(0);
        }
      }

      if (dims.size() == 0)
        dims.resize(n_dims);

      for (size_t i = 0; i < n_dims; i++)
        dims[i] = boost::lexical_cast<double>(line_pieces[i]);

      Point<T> point(dims);
      points_.push_back(dims);
      ++n_points_;
    }
  }

  ifs.close();
  std::cout << "[INFO] " << n_points_ << " points loaded with " << n_dims <<
    " dimensions each." << std::endl;
}

template<typename T>
Point<T>& Table<T>::operator[](const size_t i)
{
  if (i >= n_points_) 
  {
    std::cout << "[ERROR] Table index out of range" << std::endl;
    exit(1);
  }
  return points_[i];
}

template<typename T>
const Point<T>& Table<T>::operator[](const size_t i) const
{
  if (i >= n_points_) 
  {
    std::cout << "[ERROR] Table index out of range" << std::endl;
    exit(1);
  }
  return points_[i];
}

template<typename T>
Table<T>& Table<T>::operator=(const Table<T>& table) 
{
  n_points_ = table.n_points_;
  points_ = table.points_;
}

}; // namespace
#endif
