/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TAMatrix.hpp
  \class TAMatrix<T>
  \brief A template class for general matrices, including data storage and
  various kinds of matrix operation. This is the deinition file for the member
  methods
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/09
  \date Last modified: 2020/10/04, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <cmath>
#include "TAException.h"

using std::ios_base;
using std::cout;
using std::endl;
using std::setw;
using std::move;

template<class T>
bool inline isBasic();

template<class T>
TAMatrix<T>::TAMatrix() // the default constructor
 : fData(nullptr), fNRow(0), fNColumn(0){
} // end of the constructor

// data[i,j] = data[i*ncols+j]
template<class T>
TAMatrix<T>::TAMatrix(int nrows, int ncols, const T *data){
  if(!nrows || !ncols){
		TAMatrix();
		return;
	}
  fNRow = nrows; fNColumn = ncols;
  fData = new T[fNRow*fNColumn];
  fRowVEC.resize(fNRow); fColVEC.resize(fNColumn);

  vec_t<T> vTmp; vTmp.resize(fNColumn); // used as row vectors
  for(int i = 0; i < fNRow; i++){
		for(int j = 0; j < fNColumn; j++){
			if(data) fData[i*fNColumn+j] = data[i*fNColumn+j];
      vTmp.at(j) = fData+i*fNColumn+j;
		} // end for over columns
    fRowVEC[i] = new vec_t<T>(vTmp);
	} // end for over rows
  if(!data) Initialize();
  // assign fColVEC //
  vTmp.resize(fNRow); // used as row vectors
  for(int i = 0; i < fNColumn; i++){
    for(int j = 0; j < fNRow; j++) vTmp.at(j) = fData+j*fNColumn+i;
    fColVEC[i] = new vec_t<T>(vTmp);
	} // end for over columns
} // end of the constructor

// the copy costructor
// call default ctor here to set fData nullptr so as to
// protect it from purging in the assignment constructor
template<class T>
TAMatrix<T>::TAMatrix(const TAMatrix<T> &ma) : TAMatrix(){
  *this = ma;
} // end of the copy constructor

// move costructor
template<class T>
TAMatrix<T>::TAMatrix(TAMatrix<T> &&ma){
  fNRow = ma.nr(); fNColumn = ma.nc();
  fData = ma.fData; ma.fData = nullptr;
  fRowVEC = ma.fRowVEC; fColVEC = ma.fColVEC;
} // end of the copy constructor

template<class T>
TAMatrix<T>::~TAMatrix(){
  if(fData){ delete [] fData; fData = nullptr; }
  for(vec_t<T> *&p : fRowVEC) if(p){ delete p; p = nullptr; }
  fRowVEC.clear();
  for(vec_t<T> *&p : fColVEC) if(p){ delete p; p = nullptr; }
  fColVEC.clear();
} // end of the destructor

template<class T>
void TAMatrix<T>::Resize(int nrow, int ncol){
  if(fNRow == nrow && fNColumn == ncol) return;

  fNRow = nrow; fNColumn = ncol;
  if(fData){ delete [] fData; fData = nullptr; }
  fData = new T[fNRow*fNColumn];
  for(vec_t<T> *&p : fRowVEC) if(p){ delete p; p = nullptr; }
  for(vec_t<T> *&p : fColVEC) if(p){ delete p; p = nullptr; }
  fRowVEC.resize(fNRow);
  fColVEC.resize(fNColumn);

  vec_t<T> vTmp; vTmp.resize(fNColumn); // used as row vectors
  for(int i = 0; i < fNRow; i++){
    for(int j = 0; j < fNColumn; j++) vTmp.at(j) = fData+i*fNColumn+j;
    fRowVEC[i] = new vec_t<T>(vTmp);
  } // end for over rows
  // assign fColVEC //
  vTmp.resize(fNRow); // used as row vectors
  for(int i = 0; i < fNColumn; i++){
    for(int j = 0; j < fNRow; j++) vTmp.at(j) = fData+j*fNColumn+i;
    fColVEC[i] = new vec_t<T>(vTmp);
  } // end for over columns
} // end of member function Resize

// assignment constructor
template<class T>
TAMatrix<T> &TAMatrix<T>::operator=(const TAMatrix<T> &ma){
  if(&ma == this) return *this;

  Resize(ma.nr(), ma.nc());
  for(int i = 0; i < fNRow; i++){
		for(int j = 0; j < fNColumn; j++) (*fRowVEC[i])[j] = ma[i][j];
	} // end for over rows

  return *this;
} // end the assignment constructor

// move assignment constructor
template<class T>
TAMatrix<T> &TAMatrix<T>::operator=(TAMatrix<T> &&ma){
  fNRow = ma.nr(); fNColumn = ma.nc();
  if(fData) delete [] fData;
  fData = ma.fData; ma.fData = nullptr;

  for(vec_t<T> *&p : fRowVEC) if(p){ delete p; p = nullptr; }
  for(vec_t<T> *&p : fColVEC) if(p){ delete p; p = nullptr; }
  fRowVEC = move(ma.fRowVEC);
  fColVEC = move(ma.fColVEC);

  return *this;
} // end of move assignment constructor

// {} initialization
template<class T>
TAMatrix<T> &TAMatrix<T>::operator=(const std::initializer_list<T> &li){
  const int n = fNColumn * fNRow;
  if(int(li.size()) != n){
    TAException::Error("TAMatrix<T>", "operator={}: initializer_list size\
doesn't match the number of the matrix element.");
  }

  int cnt = 0;
  for(auto t : li){
    fData[cnt++] = t;
    if(cnt >= n) break;
  } // end for
  return *this;
} // end of the assignment constructor by an initializer list

// initilize to the unit matrix * val
template<class T>
TAMatrix<T> &TAMatrix<T>::operator=(double b){
  if(IsEmpty())
    TAException::Error("TAMatrix<T>", "operator=(double): The matrix is empty.");
  if(fNRow != fNColumn)
    TAException::Error("TAMatrix<T>", "operator=(double): The matrix is not square.");

  for(int i = 0; i < fNRow; i++){
    for(int j = 0; j < fNColumn; j++){
      if(i == j) (*fRowVEC[i])[j] = b;
      else  (*fRowVEC[i])[j] = 0;
    } // end loop over columns
  } // end of loop over rows

  return *this;
} // end of member function operator=(double)

template<class T>
vec_t<T> &TAMatrix<T>::operator[](int r){
  if(r < 0 || r >= fNRow)
    TAException::Error("TAMatrix<T>",
      "operator[]: Input row %d out of range, max: %d", r, fRowVEC.size()-1);
  return *fRowVEC[r];
} // end of member function operator[]
template<class T>
const vec_t<T> &TAMatrix<T>::operator[](int r) const{
  // XXX: return (*this)[row]; WRONG: const function won't call non-const ones
    TAException::Error("TAMatrix<T>",
      "operator[]: Input row %d out of range, max: %d", r, fRowVEC.size()-1);
  return *fRowVEC[r];
} // end of member function operator[]
template<class T>
vec_t<T> &TAMatrix<T>::cv(int c){
  if(c < 0 || c >= fNColumn)
    TAException::Error("TAMatrix<T>",
      "cv: Input column %d out of range, max: %d", c, fColVEC.size()-1);
  return *fColVEC[c];
} // end of member function operator[]
template<class T>
const vec_t<T> &TAMatrix<T>::cv(int c) const{
  // XXX: return (*this)[row]; WRONG: const function won't call non-const ones
  if(c < 0 || c >= fNColumn)
    TAException::Error("TAMatrix<T>",
      "cv: Input column %d out of range, max: %d", c, fColVEC.size()-1);
  return *fColVEC[c];
} // end of member function operator[]

// operations //
/// calculations in place
template<class T>
TAMatrix<T> &TAMatrix<T>::operator+=(const TAMatrix<T> &ma){
  Add(ma, *this); return *this;
}
template<class T>
TAMatrix<T> &TAMatrix<T>::operator-=(const TAMatrix<T> &ma){
  Subtract(ma, *this); return *this;
}
template<class T>
TAMatrix<T> &TAMatrix<T>::operator-=(const T &b){
  if(!IsSquare()) TAException::Error("TAMatrix<T>", "operator-=: The matrix is not square.");

  for(int i = 0; i < fNRow; i++) (*this)[i][i] -= b;
  return *this;
} // end of the member function operator-=
template<class T>
TAMatrix<T> TAMatrix<T>::operator-(const T &b) const{
  TAMatrix<T> ma_t((*this)); return ma_t -= b;
} // end of the member function operator-=

/// \retval returns fRowVEC * ma, NOT ma * fRowVECnian.h:23:
template<class T>
TAMatrix<T> &TAMatrix<T>::operator*=(const TAMatrix<T> &ma){
  DotProduct(ma, *this); return *this;
} // end of member function operator*

template<class T>
TAMatrix<T> &TAMatrix<T>::operator/=(const T &b){
  for(vec_t<T> *&v : fRowVEC) *v /= b;
  return *this;
} // end of member function operator/=

template<class T>
TAMatrix<T> &TAMatrix<T>::operator*=(const T &b){
  for(vec_t<T> *&v : fRowVEC) *v *= b;
  return *this;
} // end of member function operator*=

/// calculations not in place, new memory would be claimed
template<class T>
TAMatrix<T> TAMatrix<T>::operator+(const TAMatrix<T> &ma) const{
  TAMatrix<T> ma_t = *this; return ma_t += ma;
} // end operator+()

template<class T>
TAMatrix<T> TAMatrix<T>::operator-(const TAMatrix<T> &ma) const{
  TAMatrix<T> ma_t = *this; return ma_t -= ma;
} // end operator-()

template<class T>
TAMatrix<T> TAMatrix<T>::operator*(const TAMatrix<T> &ma) const{
  TAMatrix<T> ma_t(fNRow, ma.nc());
  DotProduct(ma, ma_t); return ma_t;
} // end operator*(const TAMatrix<T> &)
template<class T>
void TAMatrix<T>::DotProduct(const TAMatrix<T> &ma, TAMatrix<T> &r){ ///<\retval r=this*ma
  if(fNColumn != ma.nr())
    TAException::Error("TAMatrix<T>", "DotProduct: Dimension mismatch: ma v.s. this");
  if(!ma.fData) TAException::Error("TAMatrix<T>", "DotProduct: Input matrix fData is NULL.");

  const int ir = nr(), ic = ma.nc(), n = ma.nr();
  if(r.nr() != ir || r.nc() != ic) TAException::Error("TAMatrix<T>", "DotProduct: r's Dimension Match");
  for(int i = 0; i < ir; i++) for(int j = 0; j < ic; j++){
    r[i][j] = 0;
    for(int k = 0; k < n; k++) r[i][j] += (*this)[i][k] * ma[k][j];
  } // end for over rows and columns
} // end of member function DotProduct
template<class T> /// \retval r=this+k*ma
void TAMatrix<T>::Add(const TAMatrix<T> &ma, TAMatrix<T> &r){
  if(!ma.fData) TAException::Error("TAMatrix<T>", "Add: Input matrix fData is NULL.");
  if(DimensionMatch(ma))
    TAException::Error("TAMatrix<T>", "Add: Dimension Match");

  for(int i = 0; i < fNRow; i++) r[i] += ma[i]; // vector addition
} // end of member function Add
template<class T>
void TAMatrix<T>::Subtract(const TAMatrix<T> &ma, TAMatrix<T> &r){ ///<\retval r=this-ma
  if(!ma.fData) TAException::Error("TAMatrix<T>", "Subtract: Input matrix fData is NULL.");
  if(DimensionMatch(ma))
    TAException::Error("TAMatrix<T>", "Subtract: Dimension Match");

  for(int i = 0; i < fNRow; i++) r[i] -= ma[i]; // vector subtraction
} // end of member function Add

template<class T>
TAMatrix<T> TAMatrix<T>::operator*(const T &val) const{
  TAMatrix<T> ma_t = *this;
  return ma_t *= val;
} // end operator*(const T &)

template<class T>
TAMatrix<T> operator*(const T &b, const TAMatrix<T> &ma){
  TAMatrix<T> ma_t(ma.nr(), ma.nc());
  for(int i = 0; i < ma.nr(); i++) ma_t[i] = b*ma[i];
  return ma_t;
} // end of template function operator*

template<class T>
TAMatrix<T>::operator T() const{
  if(1 != fNRow || 1 != fNColumn)
    TAException::Error("TAMatrix<T>", "operator T: matrix not of 1x1 form.");
  return (*fRowVEC[0])[0];
} // end of member function operator <T>

// NOT inplace
template<class T>
TAMatrix<T> TAMatrix<T>::Transpose() const{
  TAMatrix<T> ma_t(fNColumn, fNRow);
  for(int i = 0; i < fNColumn; i++){
    for(int j = 0; j < fNRow; j++){
      ma_t[i][j] = (*fColVEC[i])[j]; // (*fRowVEC[j])[i] = (*fColVEC[i])[j];
    } // end for over columns
  } // end for over rows
  return ma_t;
} // end of member function Transpose

// set all the elements to zero
template<class T>
void TAMatrix<T>::Initialize(){
  if(!isBasic<T>()) return;
  for(vec_t<T> *&v : fRowVEC) for(T *t : *v) *t = 0;
} // end of member function Initialize

// display the matrix in matrix form
template<class T>
void TAMatrix<T>::Print() const{
  if(!fNRow || !fNColumn || !fData){
    TAException::Warn("TAMatrix<T>", "Matrix is empty.");
	}
	ios_base::fmtflags initial;
	initial = cout.setf(ios_base::fixed, ios_base::floatfield);
	cout.precision(6); // 15 3
	cout.unsetf(ios_base::fixed);
  cout << std::right;
	cout << "column"; // 6
	for(int i = 0; i < fNColumn; i++)
		cout << setw(10) << i;
	cout << endl;
	for(int i = 0; i < fNRow; i++){
		cout << "row" << setw(3) << i;
		for(int j = 0; j < fNColumn; j++){
			cout << "\033[32;1m" << setw(10) << (*fRowVEC[i])[j] << "\033[0m";
		} // end for over columns
		cout << endl;
//		cout << "},\n{";
	} // end for over rows
	cout.setf(initial); // restore initial formatting state
  cout << endl;
//  getchar(); // DEBUG
} // end member function Print
///< display the matrix in C/C++ readable form
template<class T>
void TAMatrix<T>::PrintInC() const{
  if(!fNRow || !fNColumn || !fData){
    TAException::Warn("TAMatrix<T>", "Matrix is empty.");
	}
	ios_base::fmtflags initial;
	initial = cout.setf(ios_base::fixed, ios_base::floatfield);
	cout.precision(8); // 15 3
	cout.unsetf(ios_base::fixed);
  cout << std::right;
	cout << endl << "{";
	for(int i = 0; i < fNRow; i++){
		cout << "{";
		for(int j = 0; j < fNColumn; j++){
			cout << "\033[32;1m" << setw(10) << (*fRowVEC[i])[j];
      if(j != fNColumn - 1) cout << ",";
      cout << "\033[0m";
		} // end for over columns
		cout << "}";
    if(i != fNRow - 1) cout << ",";
    cout << endl;
//		cout << "},\n{";
	} // end for over rows
  cout << "};" << endl;
	cout.setf(initial); // restore initial formatting state
  cout << endl;
//  getchar(); // DEBUG
}

template<class T>
bool TAMatrix<T>::IsSymmetric() const{
  if(!IsSquare()) return false;
  for(int i = 0; i < fNRow; i++){
    for(int j = 0; j < fNColumn; j++){
      if(fabs((*fRowVEC[i])[j] - (*fColVEC[i])[j]) > 1E-6) return false;
    } // end for over j
  } // end for over i
  return true;
} // end of member function IsSymmetric

/// as the name indicates. \param c must point to an object created the heap
template<class T>
void TAMatrix<T>::PushBackRow(vec_t<T> *r){
  const int n = fNColumn;
  if(n != int(r->size()))
    TAException::Error("TAMatrix<T>", "PushBackRow: Dimension mismatch.");
  fRowVEC.push_back(r); fNRow++;
  for(int i = 0; i < fNColumn; i++) fColVEC[i]->push_back(r->at(i));
} // end of member function PushBackColumn

/// as the name indicates. \param c must point to an object created the heap
template<class T>
void TAMatrix<T>::PushBackColumn(vec_t<T> *c){
  const int n = fNRow;
  if(n != int(c->size()))
    TAException::Error("TAMatrix<T>", "PushBackColumn: Dimension mismatch.");
  fColVEC.push_back(c); fNColumn++;
  for(int i = 0; i < fNRow; i++) fRowVEC[i]->push_back(c->at(i));
} // end of member function PushBackColumn
