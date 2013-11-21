/**
 * @file bmst/mlpack_code/data.hpp
 * @author Pari Ram (p.ram@gatech.edu)
 *
 * This file tests the functionalities available to the data vector 
 * and table classes
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

int main(int argc, char* argv[])
{

  std::default_random_engine generator(time(NULL));
  std::uniform_real_distribution<double> randu(0, 10);

  std::vector<double> a;
  for (size_t i = 0; i < 640000; i++)
    a.push_back(randu(generator));

  std::vector<double> b;
  for (size_t i = 0; i < 640000; i++)
    b.push_back(randu(generator));

  std::vector<double> c;
  for (size_t i = 0; i < 3000000; i++)
    c.push_back(randu(generator));

  std::vector<double> d;
  for (size_t i = 0; i < 3000000; i++)
    d.push_back(randu(generator));

  std::cout << "Testing initializer(std::vector) + operator []"; 
  const bmst::Point<double> p(a);
  const bmst::Point<double> q(b);
  const bmst::Point<double> r(c);
  const bmst::Point<double> s(d);

  assert(p.n_dims() == a.size());
  for (size_t i = 0; i < a.size(); i++) {
    assert(a[i] == p[i]);
    // std::cout << "a[i]: " << a[i] << ", p[i]: " << p[i] << std::endl;
  }

  assert(s.n_dims() == d.size());
  for (size_t i = 0; i < d.size(); i++) {
    assert(d[i] == s[i]);
    // std::cout << "d[i]: " << d[i] << ", s[i]: " << s[i] << std::endl;
  }

  std::cout << " ... PASSED" << std::endl;
  std::cout<< "Testing initializer(Point) + operator=";
  bmst::Point<double> t(p);
  bmst::Point<double> u;
  u = p;

  assert(t.n_dims() == p.n_dims());
  for (size_t i = 0; i < p.n_dims(); i++) {
    assert(p[i] == t[i]);
  }

  assert(u.n_dims() == p.n_dims());
  for (size_t i = 0; i < p.n_dims(); i++) {
    assert(p[i] == u[i]);
  }
  std::cout << " ... PASSED" << std::endl;
  std::cout << "Testing initializer() + zeros() + ones()";
  bmst::Point<double> v;
  v.zeros(10);
  assert(v.n_dims() == 10);
  for (size_t i = 0; i < 10; i++) {
    assert(v[i] == 0);
  }

  v.ones();
  assert(v.n_dims() == 10);
  for (size_t i = 0; i < 10; i++) {
    assert(v[i] == 1);
  }

  v.ones(15);
  assert(v.n_dims() == 15);
  for (size_t i = 0; i < 15; i++) {
    assert(v[i] == 1);
  }

  v.zeros();
  assert(v.n_dims() == 15);
  for (size_t i = 0; i < 15; i++) {
    assert(v[i] == 0);
  }
  std::cout << " ... PASSED" << std::endl;
  std::cout << "Testing += and -=";
  t = p;
  u = r;
  t += q;
  assert(t.n_dims() == p.n_dims());
  for (size_t i = 0; i < p.n_dims(); i++) {
    assert(t[i] == p[i] + q[i]);
  }
  u -= s;
  assert(u.n_dims() == r.n_dims());
  for (size_t i = 0; i < r.n_dims(); i++) {
    assert(u[i] == r[i] - s[i]);
  }
  std::cout << " ... PASSED" << std::endl;
  std::cout << "Testing *= and /=";
  t = p;
  u = r;
  t *= 2;
  assert(t.n_dims() == p.n_dims());
  for (size_t i = 0; i < p.n_dims(); i++) {
    assert(t[i] == p[i] * 2);
  }
  u /= 2;
  assert(u.n_dims() == r.n_dims());
  for (size_t i = 0; i < r.n_dims(); i++) {
    assert(u[i] == r[i] / 2);
  }
  std::cout << " ... PASSED" << std::endl;
  std::cout << "Testing binary + , binary -";
  t = r + s;
  u = p - q;
  assert(t.n_dims() == r.n_dims());
  for (size_t i = 0; i < r.n_dims(); i++) {
    assert(t[i] == r[i] + s[i]);
  }

  assert(u.n_dims() == p.n_dims());
  for (size_t i = 0; i < p.n_dims(); i++) {
    assert(u[i] == p[i] - q[i]);
  }
  std::cout << " ... PASSED" << std::endl;
  std::cout << "Testing binary *= , binary /=";
  t = q * 2;
  u = 5 * r;
  assert(t.n_dims() == q.n_dims());
  for (size_t i = 0; i < q.n_dims(); i++) {
    assert(t[i] == q[i] * 2);
  }

  assert(u.n_dims() == r.n_dims());
  for (size_t i = 0; i < r.n_dims(); i++) {
    assert(u[i] == r[i] * 5);
  }

  t = r / 9;
  assert(t.n_dims() == r.n_dims());
  for (size_t i = 0; i < r.n_dims(); i++) {
    assert(t[i] == r[i] / 9);
  }
  std::cout << " ... PASSED" << std::endl;
  std::cout << "Testing binary dot-product";
  double dp = 0;
  for (size_t i = 0; i < p.n_dims(); i++) 
    dp += (p[i] * q[i]);

  assert(dp == bmst::Dot(p, q));

  dp = 0;
  for (size_t i = 0; i < r.n_dims(); i++) 
    dp += (r[i] * s[i]);

  assert(dp == bmst::Dot(r, s));
  std::cout << " ... PASSED" << std::endl;


  // std::cout << "Testing failures ... ";
  // // double x = p[10];
  // // double x = r[11];
  // // t = r;
  // // t += p;
  // // t -= p;
  // // dp = bmst::Dot(r, p);

  // std::cout << " .. DONE " << std::endl;

  return 0;
}
