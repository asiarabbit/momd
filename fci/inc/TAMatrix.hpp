/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TAMatrix.hpp
  \class TAMatrix<T>
  \brief A template class for general matrices, including data storage and
  various kinds of matrix operation. This is the deinition file for the member
  methods
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/09
  \date Last modified: 2020/10/07, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iomanip>
#include <typeinfo>
#include <cmath>
#include <utility>
#include "TAException.h"

using std::ios_base;
using std::cout;
using std::endl;
using std::setw;
using std::move;

template <class T>
bool inline isBasic(); // is basic type or not, e.g., not struct or class object
template <class T>
TAMatrix<T>::TAMatrix() : fNRow(0), fNColumn(0){ // the default constructor
  fRowVEC.clear(); fRowVEC.reserve(4);
  fColVEC.clear(); fColVEC.reserve(4);
} // end of the constructor

// data[i,j] = data[i*ncols+j]
template <class T>
TAMatrix<T>::TAMatrix(int nrow, int ncol, const T *data){
  if(!nrow || !ncol){ TAMatrix(); return; }
  fNRow = nrow; fNColumn = ncol;
  fRowVEC.resize(nrow); fColVEC.resize(ncol);

  // fill the rows with dynamically allocated memories //
  for(int i = 0; i < nrow; i++){
    fRowVEC[i] = new vec_t<T>(ncol);
		if(data) for(int j = 0; j < ncol; j++) (*fRowVEC[i])[j] = data[i*ncol+j];
	} // end for over rows
  if(!data) Initialize();

  // fill the columns, but only with pointers //
  vector<T *> vc(nrow);
  for(int i = 0; i < ncol; i++){
    for(int j = 0; j < nrow; j++) vc[j] = fRowVEC[j]->at(i);
    fColVEC[i] = new vec_t<T>(vc);
	} // end for over columns
} // end of the constructor

// the copy costructor
template <class T>
TAMatrix<T>::TAMatrix(const TAMatrix<T> &ma) : TAMatrix(){ *this = ma; }

// move costructor
template <class T>
TAMatrix<T>::TAMatrix(TAMatrix<T> &&ma){ *this = ma; }
template <class T>
TAMatrix<T>::~TAMatrix<T>(){ FreeMemory(); }

template <class T>
void TAMatrix<T>::FreeMemory(){
  for(auto p : fRowVEC) for(auto &t : *p){ delete t; t = nullptr; }
  for(auto p : fColVEC) for(auto &t : *p) if(t) t = nullptr;
  for(vec_t<T> *&p : fRowVEC) if(p){ delete p; p = nullptr; }
  for(vec_t<T> *&p : fColVEC) if(p){ delete p; p = nullptr; }
  fRowVEC.clear();
  fColVEC.clear();
} // end of the destructor

template <class T>
void TAMatrix<T>::Resize(int nrow, int ncol){
  if(fNRow == nrow && fNColumn == ncol) return;

  FreeMemory();
  fNRow = nrow; fNColumn = ncol;
  fRowVEC.resize(nrow);
  fColVEC.resize(ncol);

  // fill the rows with dynamically allocated memories //
  for(int i = 0; i < nrow; i++) fRowVEC[i] = new vec_t<T>(ncol);
  Initialize();
  // fill the columns, but only with pointers //
  vector<T *> vc(nrow);
  for(int i = 0; i < ncol; i++){
    for(int j = 0; j < nrow; j++) vc[j] = fRowVEC[j]->at(i);
    // pointers shared, be aware DOUBLE FREE risk upon destroying vc
    fColVEC[i] = new vec_t<T>(vc);
	} // end for over columns
} // end of member function Resize

// assignment constructor
template <class T>
TAMatrix<T> &TAMatrix<T>::operator=(const TAMatrix<T> &ma){
  if(&ma == this) return *this;

  Resize(ma.nr(), ma.nc());
  for(int i = 0; i < fNRow; i++)
		for(int j = 0; j < fNColumn; j++) (*fRowVEC[i])[j] = ma[i][j];

  return *this;
} // end the assignment constructor

// move assignment constructor
template <class T>
TAMatrix<T> &TAMatrix<T>::operator=(TAMatrix<T> &&ma){
  if(&ma == this) return *this;

  FreeMemory();
  fNRow = ma.nr(); fNColumn = ma.nc();
  fRowVEC = ma.fRowVEC; fColVEC = ma.fColVEC;
  // so that the destruction of ma won't wreck a havoc
  for(vec_t<T> *&v : ma.fRowVEC) for(T *&t : *v) t = nullptr;
  for(vec_t<T> *&v : ma.fColVEC) for(T *&t : *v) t = nullptr;
  for(vec_t<T> *&v : ma.fRowVEC) v = nullptr;
  for(vec_t<T> *&v : ma.fColVEC) v = nullptr;

  return *this;
} // end of move assignment constructor

// {} initialization
template <class T>
TAMatrix<T> &TAMatrix<T>::operator=(const std::initializer_list<T> &li){
  const int n = fNColumn * fNRow;
  if(int(li.size()) != n){
    TAException::Error("TAMatrix<T>", "operator={}: initializer_list size\
doesn't match the number of the matrix element.");
  }

  int cnt = 0;
  for(const T &t : li){
    (*this)[cnt/fNColumn][cnt%fNColumn] = t;
    if(++cnt >= n) break;
  } // end for
  return *this;
} // end of the assignment constructor by an initializer list
// initialize to E*val
template <class T>
TAMatrix<T> &TAMatrix<T>::operator=(double b){
  if(IsEmpty())
    TAException::Error("TAMatrix<T>", "operator=(double): The matrix is empty.");
  if(fNRow != fNColumn)
    TAException::Error("TAMatrix<T>", "operator=(double): The matrix is not square.");

  for(int i = nr(); i--;) for(int j = nc(); j--;)
      if(i == j) (*this)[i][j] = b; else (*this)[i][j] = 0.;
  return *this;
} // end of member function operator=(double)


template <class T>
vec_t<T> &TAMatrix<T>::operator[](int r){
  if(r < 0 || r >= fNRow)
    TAException::Error("TAMatrix<T>",
      "operator[]: Input row %d out of range, max: %d", r, fRowVEC.size()-1);
  return *fRowVEC[r];
} // end of member function operator[]
template <class T>
const vec_t<T> &TAMatrix<T>::operator[](int r) const{
  // XXX: return (*this)[row]; WRONG: const function won't call non-const ones
  if(r < 0 || r >= fNRow)
    TAException::Error("TAMatrix<T>",
      "operator[]: Input row %d out of range, max: %d", r, fRowVEC.size()-1);
  return *fRowVEC[r];
} // end of member function operator[]
template <class T>
vec_t<T> &TAMatrix<T>::cv(int c){
  if(c < 0 || c >= fNColumn)
    TAException::Error("TAMatrix<T>",
      "cv: Input column %d out of range, max: %d", c, fColVEC.size()-1);
  return *fColVEC[c];
} // end of member function operator[]
template <class T>
const vec_t<T> &TAMatrix<T>::cv(int c) const{
  // XXX: return (*this)[row]; WRONG: const function won't call non-const ones
  if(c < 0 || c >= fNColumn)
    TAException::Error("TAMatrix<T>",
      "cv: Input column %d out of range, max: %d", c, fColVEC.size()-1);
  return *fColVEC[c];
} // end of member function operator[]

// operations //
/// calculations in place
template <class T>
TAMatrix<T> &TAMatrix<T>::operator+=(const TAMatrix<T> &ma){
  if(!DimensionMatch(ma)) TAException::Error("TAMatrix<T>", "operator+=: Dimension misatch");
    for(int i = nr(); i--;) for(int j = nc(); j--;) (*this)[i][j] += ma[i][j];
  return *this;
} // end of member function operator+=
template <class T> // this += v*ma
TAMatrix<T> &TAMatrix<T>::SelfAdd(const T &v, const TAMatrix<T> &ma){
  if(!DimensionMatch(ma)) TAException::Error("TAMatrix<T>", "operator+=: Dimension misatch");
  if(!isBasic<T>()) TAException::Error("TAMatrix<T>", "operator+=: v not basic type");
  for(int i = nr(); i--;) for(int j = nc(); j--;) (*this)[i][j] += v*ma[i][j];
  return *this;
} // end of member function SelfAdd

template <class T>
TAMatrix<T> &TAMatrix<T>::operator-=(const TAMatrix<T> &ma){
  if(!DimensionMatch(ma))
    TAException::Error("TAMatrix<T>", "operator-=: Dimension misatch");
  for(int i = nr(); i--;) for(int j = nc(); j--;) (*this)[i][j] -= ma[i][j];
  return *this;
} // end member function operator-=
template <class T>
TAMatrix<T> &TAMatrix<T>::operator-=(const T &b){
  if(!IsSquare())
    TAException::Error("TAMatrix<T>", "operator-=: The matrix is not square.");
  for(int i = 0; i < fNRow; i++) (*this)[i][i] -= b;
  return *this;
} // end of the member function operator-=

/// \retval returns fRowVEC * ma, NOT ma * fRowVECnian.h:23:
template <class T>
TAMatrix<T> &TAMatrix<T>::operator*=(const TAMatrix<T> &ma){
  TAMatrix<T> prod;
  DotProduct(ma, prod); return *this = prod;
} // end of member function operator*=

template <class T>
TAMatrix<T> &TAMatrix<T>::operator/=(const T &b){
  for(vec_t<T> *&v : fRowVEC) *v /= b;
  return *this;
} // end of member function operator/=

template <class T>
TAMatrix<T> &TAMatrix<T>::operator*=(const T &b){
  for(vec_t<T> *&v : fRowVEC) *v *= b;
  return *this;
} // end of member function operator*=

template <class T>
void TAMatrix<T>::DotProduct(const TAMatrix<T> &ma, TAMatrix<T> &r) const{ ///<\retval r=this*ma
  if(fNColumn != ma.nr())
    TAException::Error("TAMatrix<T>", "DotProduct: Dimension mismatch: ma v.s. this");
  if(&ma == &r)
    TAException::Error("TAMatrix<T>", "DotProduct: the right matrix cannot be the same as the result matrix");

  const int ir = nr(), ic = ma.nc(), n = ma.nr();
  if(r.nr() != ir || r.nc() < ic)
    TAException::Error("TAMatrix<T>", "DotProduct: the result matrix is not of right dimensions");
  for(int i = 0; i < ir; i++) for(int j = 0; j < ic; j++){
    r[i][j] = 0;
    for(int k = 0; k < n; k++) r[i][j] += (*this)[i][k] * ma[k][j];
  } // end for over rows and columns
} // end of member function DotProduct
template <class T> /// \retval r=this+k*ma
void TAMatrix<T>::Add(const TAMatrix<T> &ma, TAMatrix<T> &r) const{
  if(!DimensionMatch(ma) || !DimensionMatch(r))
    TAException::Error("TAMatrix<T>", "Add: Dimension misatch");
  for(int i = nr(); i--;) for(int j = nc(); j--;) r[i][j] = (*this)[i][j] + ma[i][j];
} // end of member function Add
template <class T>
void TAMatrix<T>::Subtract(const TAMatrix<T> &ma, TAMatrix<T> &r) const{ ///<\retval r=this-ma
  if(!DimensionMatch(ma) || !DimensionMatch(r))
    TAException::Error("TAMatrix<T>", "Subtract: Dimension misatch");
  for(int i = nr(); i--;) for(int j = nc(); j--;) r[i][j] = (*this)[i][j] - ma[i][j];
} // end of member function Subtract
template <class T> // r=this*v
void TAMatrix<T>::Scale(const T &v, TAMatrix<T> &r){
  if(!DimensionMatch(r)) TAException::Error("TAMatrix<T>", "Scale: Dimension misatch");
  if(!isBasic<T>()) TAException::Error("TAMatrix<T>", "Scale: v not basic type");
  for(int i = nr(); i--;) for(int j = nc(); j--;) r[i][j] = v*(*this)[i][j];
} // end of member function Scale

template <class T>
TAMatrix<T>::operator T() const{
  if(1 != fNRow || 1 != fNColumn)
    TAException::Error("TAMatrix<T>", "operator T: matrix not of 1x1 form.");
  return (*this)[0][0];
} // end of member function operator <T>

// NOT inplace
template <class T>
TAMatrix<T> TAMatrix<T>::Transpose() const{
  TAMatrix<T> ma_t(fNColumn, fNRow);
  for(int i = nc(); i--;) for(int j = nr(); j--;) ma_t[i][j] = (*this)[j][i];
  return ma_t;
} // end of member function Transpose

// set all the elements to zero
template <class T>
void TAMatrix<T>::Initialize(){ for(vec_t<T> *v : fRowVEC) v->initialize(); }
template <class T> /// initialize to diagonal matrix with array d
void TAMatrix<T>::diag(T *d){
  if(!IsSquare()) TAException::Warn("TAMatrix<T>", "diag: Not a square matrix");
  if(!isBasic<T>()) TAException::Warn("TAMatrix<T>", "diag: T not of basic type");
  for(int i = nr(); i--;) for(int j = nc(); j--;)
    if(i == j) (*this)[i][j] = d[i];
    else (*this)[i][j] = 0.;
} // end of member function diag

/// copy the first (nr_,nc_) block
void TAMatrix<T>::copy(const TAMatrix<T> &ma, int nr_, int nc_){
  if(ma.nr() < nr_ || ma.nc() < nc_)
    TAException::Error("TAMatrix<T>", "copy: Input matrix smaller than required");
  for(int i = 0; i < nr_; i++) for(int j = 0; j < nc_; j++) this[i][j] = ma[i][j];
} // end of member function copy

// display the matrix in matrix form
template <class T>
ostream &TAMatrix<T>::Print(ostream &os) const{
  if(!fNRow || !fNColumn){
    TAException::Warn("TAMatrix<T>", "Matrix is empty.");
	}
	ios_base::fmtflags initial;
	initial = os.setf(ios_base::fixed, ios_base::floatfield);
	os.precision(6); // 15 3
	os.unsetf(ios_base::fixed);
  os << std::right;
	os << "column"; // 6
	for(int i = 0; i < fNColumn; i++)
		os << setw(15) << i;
	os << endl;
	for(int i = 0; i < fNRow; i++){
		os << "row" << setw(3) << i;
		for(int j = 0; j < fNColumn; j++){
			os << "\033[32;1m" << setw(15) << (*fRowVEC[i])[j] << "\033[0m";
		} // end for over columns
		os << endl;
//		os << "},\n{";
	} // end for over rows
	os.setf(initial); // restore initial formatting state
  return os << endl;
//  getchar(); // DEBUG
} // end member function Print
template <class T>
void TAMatrix<T>::Print() const{ Print(cout); } ///< display the matrix in matrix form

///< display the matrix in C/C++ readable form
template <class T>
void TAMatrix<T>::PrintInC() const{
  if(!fNRow || !fNColumn){
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
} // end of member function PrintInC

template <class T>
bool TAMatrix<T>::IsSymmetric() const{
  if(!IsSquare()) return false;
  for(int i = 0; i < fNRow; i++){
    for(int j = i+1; j < fNColumn; j++){
      const T g = 100.*fabs((*fRowVEC[i])[j] - (*fRowVEC[j])[i]);
      const T &d1 = (*fRowVEC[i])[i], &d2 = (*fRowVEC[j])[j];
      // less than (machine precision+2)*(corresponding diagonal elements)
      if(fabs(d1)+g != fabs(d1) || fabs(d2)+g != fabs(d2)) return false;
    } // end for over j
  } // end for over i
  return true;
} // end of member function IsSymmetric

/// as the name indicates. \param c must point to an object created the heap
template <class T>
void TAMatrix<T>::PushBackRow(const vec_t<T> &r){
  if(fNColumn != int(r.size()))
    TAException::Error("TAMatrix<T>", "PushBackRow: Dimension mismatch.");
  vec_t<T> *rv = new vec_t<T>(r);
  fRowVEC.push_back(rv); fNRow++;
  for(int i = 0; i < fNColumn; i++) fColVEC[i]->push_back(rv->at(i));
} // end of member function PushBackColumn

/// as the name indicates. \param c must point to an object created the heap
template <class T>
void TAMatrix<T>::PushBackColumn(const vec_t<T> &c){
  if(fNRow != int(c.size()))
    TAException::Error("TAMatrix<T>", "PushBackColumn: Dimension mismatch.");
  vec_t<T> *cv = new vec_t<T>(c);
  fColVEC.push_back(cv); fNColumn++;
  for(int i = 0; i < fNRow; i++) fRowVEC[i]->push_back(cv->at(i));
} // end of member function PushBackColumn

template <class T>
ostream& operator<<(ostream &os, const TAMatrix<T> &m){ return m.Print(os); }

/// (*this) = Q^T*Q*(*this), so this is more orthogonal to vectors in Q
/// this is also dubbed ** full reorthogonalization **, operation count: n*m^2
/// i.e. scales with the ietartion count n
template <class T>
void TAMatrix<T>::Purify(const TAMatrix<T> &q){
  if(q.nr() != nr()) TAException::Error("TAMatrix<T>", "Purify: Dimension mismatch.");
  const int n = q.nc(), m = nr(); T s[m]{};
  for(int k = m; k--;) for(int i = n; i--;) for(int l = m; l--;)
    s[k] += q[k][i]*q[l][i]*(*this)[l][0];
  for(int k = 0; k < m; k++) (*this)[k][0] -= s[k];
} // end of member function Purify
