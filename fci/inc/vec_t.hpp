/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file vec_t<T>.hpp
  \class vec_t<T>
  \brief Nested class of TAMatrix<T>, to represent the reference of the column
  or the row vectors of a TAMatrix<T> object.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/26
  \date Last modified: 2020/02/29, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "TAException.h"

using std::cout;
using std::endl;
using std::setw;
using std::ios_base;

template<class T>
bool inline isBasic(){
  if(typeid(T) == typeid(double) || typeid(T) == typeid(int)
      || typeid(T) == typeid(short) || typeid(T) == typeid(unsigned) ){
    return true;
  }
  return false;
} // end inline function isBasic

// vector struct for the constituent vectors in the matrix
template<class T>
vec_t<T>::vec_t(int n) : fData(nullptr){
  if(0 != n){
    fData = new T[n]; this->resize(n);
    for(int i = 0; i < n; i++) this->at(i) = fData + i;
    if(isBasic<T>()) for(T *&t : (*this)) *t = 0;
  }
}
template<class T>
vec_t<T>::vec_t(const vec_t<T> &v) : fData(nullptr){
  // construction and assignment are essentially different operations
  // for this class. Upon construction, the pointers in this vector are set
  // and fixed against assignment particularly, i.e., assignment CANNOT change
  // the pointers in **this** vector. Well, copy ctor is supposed to do so
  const int n = v.size(); this->resize(n);
  for(int i = 0; i < n; i++) this->at(i) = v.at(i);
}
template<class T>
vec_t<T>::vec_t(vec_t &&v) : fData(nullptr){
  (*this) = move(v);
}
template<class T>
vec_t<T>::~vec_t(){
  if(fData){
    delete [] fData; fData = nullptr;
  }
  this->clear();
}
// only pass values, won't change length and pointers
template<class T>
vec_t<T> &vec_t<T>::operator=(const vec_t<T> &v){
  if(v.size() != this->size()){
    TAException::Error("vec_t<T>", "operator=: Dimension mismatch.");
  }
  for(int i = v.size(); i--;) (*this)[i] = v[i];
  return *this;
}
template<class T>
vec_t<T> &vec_t<T>::operator=(vec_t &&v){
  const int n = v.size(); this->resize(n);
  for(int i = 0; i < n; i++) this->at(i) = v.at(i);

  // this is an independent vector, which owns its own data //
  if(v.fData){
    fData = v.fData;
    v.fData = nullptr;
  }

  return *this;
}
template<class T>
void vec_t<T>::SetUniformValue(const T &b){
  for(T *t : *this) *t = b;
}
template<class T>
T &vec_t<T>::operator[](int i){
  if(i < 0 || i >= int(this->size())){
    TAException::Error("vec_t<T>", "operator[]: Input i=%d, out of range.", i);
  }
  return *this->at(i);
}
template<class T>
const T &vec_t<T>::operator[](int i) const{
  /// XXX: return (*this)[i]; WRONG: trigger self-calling, an endless recursion
  if(i < 0 || i >= int(this->size())){
    TAException::Error("vec_t<T>",
      "operator[] const: Input i=%d, out of range.", i);
  }
  return *this->at(i);
}
template<class T>
vec_t<T> vec_t<T>::operator+(
    const vec_t<T> &v) const{
  const int n = this->size();
  if(v.size() != n)
    TAException::Error("vec_t<T>", "operator+: Dimension mismatch.");
  vec_t<T> sum(n);void SetUniformValue(const T &b);
  for(int i = n; i--;){
    sum[i] = (*this)[i] + v[i];
  }
  return sum;
}
template<class T>
vec_t<T> vec_t<T>::operator-(const vec_t<T> &v) const{
  const int n = this->size(), nn = v.size();
  if(n != nn) TAException::Error("vec_t<T>", "operator+: Dimension mismatch.");
  vec_t<T> sum(n);
  for(int i = 0; i < n; i++) sum[i] = (*this)[i] - v[i];
  return sum;
}
template<class T>
vec_t<T> vec_t<T>::operator*(const vec_t<T> &v) const{
  const int n = this->size(), nn = v.size();
  if(!nn || !n)
    TAException::Error("vec_t<T>", "operator*: empty vector(s) found.");

  if(nn == 1 || n == 1){ // number * vector or vector * number
    const int l = n > nn ? n : nn;
    vec_t<T> prod(l);
    for(int i = 0; i < l; i++){
      if(n > 1) prod[i] = (*this)[i] * v[0];
      else prod[i] = (*this)[0] * v[i];
    } // end for
    return prod;
  } // end outer if
  if(n != nn) TAException::Error("vec_t<T>", "operator*: Dimension mismatch.");
  bool basicType = isBasic<T>();
  vec_t<T> prod(1); // automatically initialized to zero upon construction
  for(int i = n; i--;){
    if(basicType && (!v[i] || !(*this)[i]) ) continue;
    if(1. == v[i]) prod[0] += (*this)[i];
    else if(1. == (*this)[i]) prod[0] += v[i];
    else prod[0] += (*this)[i] * v[i];
  }
  return prod;
}

template<class T>
vec_t<T> vec_t<T>::operator/(const vec_t<T> &v) const{
  if(1 != int(v.size())){
    TAException::Error("vec_t<T>", "operator/: Input vector not of 1x1 form.");
  }
  return (*this)*v[0];
}
template<class T>
vec_t<T> vec_t<T>::operator*(const T &v) const{
  const int n = this->size();
  vec_t<T> prod(n);
  for(int i = 0; i < n; i++) prod[i] = (*this)[i] * v;
  return prod;
}
template<class T>
vec_t<T> vec_t<T>::operator/(const T &v) const{
  const int n = this->size();
  vec_t<T> prod(n);
  if(isBasic<T>() && v == 0)
    TAException::Error("vec_t<T>", "operator/: Input object is zero.");
  for(int i = 0; i < n; i++) prod[i] = (*this)[i] / v;
  return prod;
}
template<class T>
vec_t<T> &vec_t<T>::operator+=(const vec_t<T> &v){
  const int n = this->size(), nn = v.size();
  if(n != nn) TAException::Error("vec_t<T>", "operator+=: Dimension mismatch.");
  for(int i = n; i--;) (*this)[i] += v[i];
  return *this;
}
template<class T>
vec_t<T> &vec_t<T>::operator-=(const vec_t<T> &v){
  const int n = this->size(), nn = v.size();
  if(n != nn) TAException::Error("vec_t<T>", "operator-=: Dimension mismatch.");
  for(int i = n; i--;) (*this)[i] -= v[i];
  return *this;
}
template<class T>
vec_t<T> &vec_t<T>::operator*=(const T &b){
  for(T *t : (*this)) *t *= b;
  return *this;
}
template<class T>
vec_t<T> &vec_t<T>::operator/=(const T &b){
  if(!isBasic<T>()){
    TAException::Error("vec_t<T>", "operator/=: Input not of basic type.");
  }
  if(!b) TAException::Error("vec_t<T>", "operator/=: Input is zero.");
  for(T *t : (*this)) *t /= b;
  return *this;
}
// change the size of the vector
template<class T>
void vec_t<T>::Resize(int n){
  if(fData) delete fData;
  fData = new T[n]; this->resize(n);
  for(int i = 0; i < n; i++) this->at(i) = fData + i;
  if(isBasic<T>()) for(T *&t : (*this)) *t = 0;
}
template<class T>
T vec_t<T>::norm2() const{
  return ((*this)*(*this))[0];
}
template<class T>
T vec_t<T>::norm() const{
  return sqrt(norm2());
}
template<class T>
vec_t<T> &vec_t<T>::normalize(){
  const T m(norm());
  if(isBasic<T>() && m == 1.) return *this;
  for(T *t : (*this)) *t /= m;
  return *this;
}
template<class T>
void vec_t<T>::Print() const{
  const int n = this->size();
  cout << "vec_t<T> Print: totally " << n << " elements." << endl;
  ios_base::fmtflags initial = cout.setf(ios_base::fixed, ios_base::floatfield);
  cout.unsetf(ios_base::floatfield);
  cout.precision(6);
  cout << std::right;
  for(int i = 0; i < n; i++) cout << setw(10) << i;
  cout << "\033[32;1m" << endl;
  for(int i = 0; i < n; i++) cout << setw(10) << (*this)[i];
  cout << "\033[0m" << endl;
  cout.setf(initial);
}
