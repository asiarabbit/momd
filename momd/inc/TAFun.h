/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFun.h
  \class TAFun
  \brief To define a functor so as to be used in template functions that require
  a functor to be passed as an argument.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/25
  \date Last modified: 2020/09/12 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAFun_h_
#define _TAFun_h_

template<class T>
class TAFun{
public:
  TAFun(){}
  ~TAFun(){}

  /// this is for general function f(x)
  virtual T operator()(double x) const{ return T(0); }
  /// this is for function vector f(n, x, f), i.e. f_i(x_0,x_1,..x_n)
  /// \param n is the dimension of x (f)
  virtual void operator()(int n, const T *x, T *f) const{}
};

#endif
