/**
 * @file bregman_mst/mlpack_code/data.hpp
 * @author Pari Ram (p.ram@gatech.edu)
 *
 * This file implements the data storage and some minor vector and 
 * matrix operations. 
 */

#ifndef BMST_DATA_HPP_
#define BMST_DATA_HPP_

#include <string>
#include <vector>

namespace bmst
{ 

template <typename T>
class Point
{
private:
  std::vector<T> values_;
  size_t n_dims_;

public:
  Point();
  Point(const std::vector<T>& point);
  Point(const Point<T>& point);

  const size_t n_dims() const { return n_dims_; }

  T& operator[](const size_t i);
  const T& operator[](const size_t i) const;
  Point<T>& operator=(const Point<T>& point);
  Point<T>& operator+=(const Point<T>& point);
  Point<T>& operator-=(const Point<T>& point);
  Point<T>& operator*=(const double scalar);
  Point<T>& operator/=(const double scalar);
  void zeros();
  void zeros(const size_t num_dims);

  void ones();
  void ones(const size_t num_dims);
  
  void print() const;
  
};


// Binary +
template<typename T>
Point<T> operator+(const Point<T>& lhs, const Point<T>& rhs);

// Binary -
template<typename T>
Point<T> operator-(const Point<T>& lhs, const Point<T>& rhs);

// Binary *
template<typename T>
Point<T> operator*(const double lhs, const Point<T>& rhs);

template<typename T>
Point<T> operator*(const Point<T>& lhs, const double rhs);

// Binary /
template<typename T>
Point<T> operator/(const Point<T>& lhs, const double rhs);

// Binary dot product
template<typename T>
double Dot(const Point<T>& a, const Point<T>& b);

template <typename T>
class Table
{
private:
  std::vector<Point<T> > points_;
  size_t n_points_;

public:
  Table();
  Table(const std::vector<std::vector<T> >& points);
  Table(const std::vector<Point<T> >& points);
  Table(const Table& table);
  // Read table in from a file
  Table(const std::string& file_name);

  const size_t n_points() const { return n_points_; }

  Point<T>& operator[](const size_t i);
  const Point<T>& operator[](const size_t i) const;
  Table& operator=(const Table& table);
  
  void print() const;
  
}; // class

}; // namespace

#include "data_impl.hpp"

#endif
