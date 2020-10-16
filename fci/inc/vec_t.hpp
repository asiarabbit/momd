/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file vec_t<T>.hpp
  \class vec_t<T>
  \brief Nested class of TAMatrix<T>, to represent the reference of the column
  or the row vectors of a TAMatrix<T> object.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/26
  \date Last modified: 2020/10/07, by SUN Yazhou
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

template <class T>
bool inline isBasic(){
  if(typeid(T) == typeid(double) || typeid(T) == typeid(int)
      || typeid(T) == typeid(float) || typeid(T) == typeid(bool)
      || typeid(T) == typeid(short) || typeid(T) == typeid(unsigned) ){
    return true;
  }
  return false;
} // end inline function isBasic

// vector struct for the constituent vectors in the matrix
template <class T>
vec_t<T>::vec_t(int n){
  if(0 != n){
    this->reserve(n);
    for(int i = n; i--;) this->push_back(new T);
    initialize();
  }
} // end of constructor
template <class T>
void vec_t<T>::initialize(){ if(isBasic<T>()) for(T *&t : *this) *t = 0.; }
template <class T>
vec_t<T>::vec_t(const vec_t &v){
  Resize(v.size());
  *this = v;
} // end of copy constructor
template <class T>
vec_t<T>::vec_t(vec_t &&v){ (*this) = move(v); }
template <class T>
void vec_t<T>::FreeMemory(){
  for(T *&p : *this) if(p){ delete p; p = nullptr; }
  this->clear();
} // end FreeMemory
template <class T>
vec_t<T>::~vec_t(){
  FreeMemory();
} // end of destructor
// only pass values, won't change length and pointers
template <class T>
vec_t<T> &vec_t<T>::operator=(const vec_t<T> &v){
  if(v.size() != this->size())
    TAException::Error("vec_t<T>", "operator=: Dimension mismatch.");
  for(int i = v.size(); i--;) (*this)[i] = v[i];
  return *this;
} // end of assignment
template <class T>
vec_t<T> &vec_t<T>::operator=(vec_t &&v){
  FreeMemory();
  const int n = v.size(); this->resize(n);
  for(int i = 0; i < n; i++) this->at(i) = v.at(i);
  for(T *&t : v) t = nullptr;
  v.clear();

  return *this;
} // end of move assignment
template <class T>
void vec_t<T>::SetUniformValue(const T &b){ for(T *t : *this) *t = b; }
template <class T>
T &vec_t<T>::operator[](int i){
  if(i < 0 || i >= int(this->size())){
    TAException::Error("vec_t<T>", "operator[]: Input i=%d, out of range.", i);
  }
  return *this->at(i);
} // end of member function operator[]
template <class T>
const T &vec_t<T>::operator[](int i) const{
  /// XXX: return (*this)[i]; WRONG: trigger self-calling, an endless recursion
  if(i < 0 || i >= int(this->size())){
    TAException::Error("vec_t<T>",
      "operator[] const: Input i=%d, out of range.", i);
  }
  return *this->at(i);
} // end of member function operator[] const
template <class T>
vec_t<T> &vec_t<T>::operator+=(const vec_t<T> &v){
  const int n = this->size(), nn = v.size();
  if(n != nn) TAException::Error("vec_t<T>", "operator+=: Dimension mismatch.");
  for(int i = n; i--;) (*this)[i] += v[i];
  return *this;
} // end of member function operator+=(const vec_t<T> &v)
template <class T>
vec_t<T> &vec_t<T>::operator-=(const vec_t<T> &v){
  const int n = this->size(), nn = v.size();
  if(n != nn) TAException::Error("vec_t<T>", "operator-=: Dimension mismatch.");
  for(int i = n; i--;) (*this)[i] -= v[i];
  return *this;
} // end of member function operator-=(const vec_t<T> &v)
template <class T>
vec_t<T> &vec_t<T>::operator*=(const T &b){
  for(T *t : (*this)) *t *= b;
  return *this;
}
template <class T>
vec_t<T> &vec_t<T>::operator/=(const T &b){
  if(!isBasic<T>()){
    TAException::Error("vec_t<T>", "operator/=: Input not of basic type.");
  }
  if(!b) TAException::Error("vec_t<T>", "operator/=: Input is zero.");
  for(T *t : (*this)) *t /= b;
  return *this;
} // end of member function operator/=(const T &b)
// change the size of the vector
template <class T>
void vec_t<T>::Resize(int n){
  if(n == int(this->size())) return;
  FreeMemory();

  this->reserve(n);
  for(int i = 0; i < n; i++) this->push_back(new T);
  if(isBasic<T>()) for(T *t : (*this)) *t = 0.;
} // end of member function Resize
template <class T>
inline T sqr(const T &t){ return t*t; }
template <class T>
T vec_t<T>::norm2() const{ T s = 0; for(T *t : (*this)) s += sqr(*t); return s; }
template <class T>
T vec_t<T>::norm() const{ return sqrt(norm2()); }
template <class T>
vec_t<T> &vec_t<T>::normalize(){
  const T m(norm());
  if(isBasic<T>() && m == 1.) return *this;
  for(T *t : (*this)) *t /= m;
  return *this;
} // end of member function normalize
template <class T>
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
} // end of member function Print
