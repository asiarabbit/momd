/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file vec_t.h
  \class TAMatrix<T>::vec_t
  \brief Nested class of TAMatrix<T>, to represent the reference of the column
  or the row vectors of a TAMatrix<T> object.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/29
  \date Last modified: 2020/02/29, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _vec_t_h_
#define _vec_t_h_

#include <vector>

using std::vector;

// vector struct for the constituent vectors in the matrix
template<class T>
struct vec_t : public vector<T *>{
  /// vector length, n != 0 is for independent vectors only
  // (not associated with a matrix)
  vec_t(int n = 0);
  vec_t(const vec_t &v);
  vec_t(vec_t &&v);
  virtual ~vec_t();

  vec_t &operator=(const vec_t &v);
  vec_t &operator=(vec_t &&v);
  void SetUniformValue(const T &b);
  T &operator[](int i);
  const T &operator[](int i) const;
  vec_t operator+(const vec_t &v) const;
  vec_t operator-(const vec_t &v) const;
  vec_t operator*(const vec_t &v) const;
  vec_t operator/(const vec_t &v) const;
  vec_t operator*(const T &v) const;
  vec_t operator/(const T &v) const;
  vec_t &operator+=(const vec_t &v);
  vec_t &operator-=(const vec_t &v);
  vec_t &operator*=(const T &b);
  vec_t &operator/=(const T &b);
  void Resize(int n); // change the size of the vector
  T norm() const;
  T norm2() const; // norm^2
  vec_t &normalize();
  void Print() const;

  template<class T1>
  friend vec_t<T1> operator*(const T1 &b, const vec_t<T1> &v);

  T *fData; ///< an expedient approach to accommodate rvalue objects
};


#include "vec_t.hpp"

#endif
