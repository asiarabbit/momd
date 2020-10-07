/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TAMatrix.h
  \class TAMatrix<T>
  \brief a template class for general matrices, including data storage and
  various kinds of matrix operations.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/09
  \date Last modified: 2020/10/07, by SUN Yazhou
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

template<class T>
class TAMatrix{
public:
  TAMatrix();
  /// data[i,j] = data[i*ncols+j]
	TAMatrix(int nrows, int ncols, const T *data = nullptr);
  explicit TAMatrix(int nrows) : TAMatrix(nrows, 1){} ///< a vector
	TAMatrix(const TAMatrix<T> &ma); ///< the copy costructor
  TAMatrix(TAMatrix<T> &&ma); // move constructor
	virtual ~TAMatrix();

  TAMatrix<T> &operator=(const TAMatrix<T> &ma); ///< assignment constructor
  TAMatrix<T> &operator=(TAMatrix<T> &&ma); ///< move assignment constructor
  /// {} initialization
  TAMatrix<T> &operator=(const std::initializer_list<T> &li);
  TAMatrix<T> &operator=(double val); ///< initialize to E*val
  /// initialize to E*val
  TAMatrix<T> &operator=(int val){ return (*this) = double(val); }
  vec_t<T> &operator[](int row); ///< operator[row][column]
  const vec_t<T> &operator[](int row) const; ///< const version
  bool DimensionMatch(const TAMatrix &b) const{ return b.nr() == nr() && b.nc() == nc(); }

  // operations //
  /// calculations in place
	TAMatrix<T> &operator+=(const TAMatrix<T> &ma);
	TAMatrix<T> &operator-=(const TAMatrix<T> &ma);
	TAMatrix<T> &operator-=(const T &b); ///< -b*I
  /// \retval returns (*this) * ma, NOT ma * (*this)
	TAMatrix<T> &operator*=(const TAMatrix<T> &ma);
  /// \retval returns (*this) * val, NOT val * (*this)
	TAMatrix<T> &operator*=(const T &val);
  TAMatrix<T> &operator/=(const T &val); // only valid for T==double
  /// calculations not in place, new memory would be claimed
	TAMatrix<T> operator+(const TAMatrix<T> &ma) const;
	TAMatrix<T> operator-(const TAMatrix<T> &ma) const;
	TAMatrix<T> operator-(const T &b) const; ///< -b*I
	TAMatrix<T> operator*(const TAMatrix<T> &ma) const;
  TAMatrix<T> operator*(const T &val) const;
  void DotProduct(const TAMatrix<T> &ma, TAMatrix<T> &r) const; ///<\retval r=this*ma
  void Add(const TAMatrix<T> &ma, TAMatrix<T> &r) const; ///<\retval r=this+ma
  void Subtract(const TAMatrix<T> &ma, TAMatrix<T> &r) const; ///<\retval r=this-ma

  operator T() const;
	TAMatrix<T> Transpose() const; ///< NOT inplace
  void Initialize(); ///< set all the elements to zero
  // re-shape the matrix, do nothing if the shape remains
  void Resize(const int nrow, const int ncol);

  int nr() const{ return fNRow; }
  int nc() const{ return fNColumn; }
  vec_t<T> &rv(int r){ return (*this)[r]; } ///< row vector i
  const vec_t<T> &rv(int r) const{ return (*this)[r]; } ///< const version
  vec_t<T> &cv(int c); ///< col vector i
  const vec_t<T> &cv(int c) const; ///< const version

  void Print() const; ///< display the matrix in matrix form
  void PrintInC() const; ///< display the matrix in C/C++ readable form
  bool IsEmpty() const{ return 0 == nr() || 0 == nc(); }
  bool IsVector() const{ return nc() == 1; }
  bool IsSquare() const{ return nr() == nc(); }
  bool IsSymmetric() const;
  /// as the name indicates. \param c(r) must point to an object in the heap
  void PushBackRow(vec_t<T> *r);
  void PushBackColumn(vec_t<T> *c);

  template<class T1>
  friend TAMatrix<T1> operator*(const T1 &val, const TAMatrix<T1> &ma);

private:
  T *fData;
  int fNRow;
  int fNColumn;
  vector<vec_t<T> *> fRowVEC;
  vector<vec_t<T> *> fColVEC;
};


typedef TAMatrix<double> matrix;
typedef TAMatrix<double> TAMatrix2D;
typedef TAMatrix<TAMatrix2D> TAMatrix4D; /// 4-th order tensor
typedef TAMatrix<TAMatrix4D> TAMatrix6D; /// 6-th order tensor

#include "TAMatrix.hpp" // the definition of template class TAMatri

#endif
