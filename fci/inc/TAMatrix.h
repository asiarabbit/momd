/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TAMatrix.h
  \class TAMatrix<T>
  \brief a template class for general matrices, including data storage and
  various kinds of matrix operations.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/09
  \date Last modified: 2020/10/08, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMatrix_h_
#define _TAMatrix_h_

#include <iostream>
#include <initializer_list>

#include "vec_t.h"

using std::vector;
using std::ostream;

template <class T>
class TAMatrix{
public:
  TAMatrix();
  /// data[i,j] = data[i*ncols+j]
	TAMatrix(int nrow, int ncol, const T *data = nullptr);
  explicit TAMatrix(int nrows) : TAMatrix(nrows, 1){} ///< a vector
	TAMatrix(const TAMatrix<T> &ma); ///< the copy costructor
  TAMatrix(TAMatrix<T> &&ma); // move constructor
	virtual ~TAMatrix();
  void FreeMemory(); // an destructor agent //

  TAMatrix<T> &operator=(const TAMatrix<T> &ma); ///< assignment constructor
  TAMatrix<T> &operator=(TAMatrix<T> &&ma); ///< move assignment constructor
  /// {} initialization
  TAMatrix<T> &operator=(const std::initializer_list<T> &li);
  TAMatrix<T> &operator=(double v); ///< initialize to E*val
  TAMatrix<T> &operator=(int v){ return (*this) = double(v); } ///< initialize to E*val
  /// initialize to E*val
  vec_t<T> &operator[](int row); ///< operator[row][column]
  const vec_t<T> &operator[](int row) const; ///< const version
  bool DimensionMatch(const TAMatrix &b) const{ return b.nr() == nr() && b.nc() == nc(); }

  // operations //
  /// calculations in place
	TAMatrix<T> &operator+=(const TAMatrix<T> &ma);
	TAMatrix<T> &SelfAdd(const T &v, const TAMatrix<T> &ma); ///< this += v*ma
	TAMatrix<T> &operator-=(const TAMatrix<T> &ma);
	TAMatrix<T> &operator-=(const T &b); ///< -b*I
  /// \retval returns (*this) * ma, NOT ma * (*this)
	TAMatrix<T> &operator*=(const TAMatrix<T> &ma);
  /// \retval returns (*this) * val, NOT val * (*this)
	TAMatrix<T> &operator*=(const T &val);
  TAMatrix<T> &operator/=(const T &val); // only valid for T==double
  void DotProduct(const TAMatrix<T> &ma, TAMatrix<T> &r) const; ///<\retval r=this*ma
  void Add(const TAMatrix<T> &ma, TAMatrix<T> &r) const; ///<\retval r=this+ma
  void Subtract(const TAMatrix<T> &ma, TAMatrix<T> &r) const; ///<\retval r=this-ma
  void Scale(const T &v, TAMatrix<T> &r); ///< r=v*this

  explicit operator T() const;
	TAMatrix<T> Transpose() const; ///< NOT inplace
  void Initialize(); ///< set all the elements to zero
  void diag(T *d); /// initialize to diagonal matrix with array d
  void copy(const TAMatrix<T> &ma, int nr_, int nc_); // copy the first (nr_,nc_) block
  // re-shape the matrix, do nothing if the shape remains
  void Resize(const int nrow, const int ncol);

  int nr() const{ return fNRow; }
  int nc() const{ return fNColumn; }
  vec_t<T> &rv(int r){ return (*this)[r]; } ///< row vector i
  const vec_t<T> &rv(int r) const{ return (*this)[r]; } ///< const version
  vec_t<T> &cv(int c); ///< col vector i
  const vec_t<T> &cv(int c) const; ///< const version

  void Print() const; ///< display the matrix in matrix form
  ostream &Print(ostream &os) const; ///< display the matrix in matrix form into os
  void PrintInC() const; ///< display the matrix in C/C++ readable form
  bool IsEmpty() const{ return 0 == nr() || 0 == nc(); }
  bool IsVector() const{ return nc() == 1; }
  bool IsSquare() const{ return nr() == nc(); }
  bool IsSymmetric() const;
  /// as the name indicates. these methods push a copy of r & c to the matrix
  void PushBackRow(const vec_t<T> &r);
  void PushBackColumn(const vec_t<T> &c);
  void EraseColumn(int ic0, int ic1); // remove column [ic0, ic1)
  /// (*this) = Q^T*Q*(*this), so this is more orthogonal to vectors in Q
  void Purify(const TAMatrix<T> &Q);

  template <class T1>
  friend ostream& operator<<(ostream &os, const TAMatrix<T1> &m);

private:
  int fNRow, fNColumn;
  vector<vec_t<T> *> fRowVEC;
  vector<vec_t<T> *> fColVEC;
};


typedef TAMatrix<double> matrix;

#include "TAMatrix.hpp" // the definition of template class TAMatrix

#endif
