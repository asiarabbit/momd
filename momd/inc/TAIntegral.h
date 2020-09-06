/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAIntegral.h
  \class TAIntegral
  \brief This class implements general numerical integral. Input and results are
  all stored in arrays. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/22 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAIntegral_h_
#define _TAIntegral_h_

#include "TAException.h"

template<typename T, typename FUNCTOR>
class TAIntegral{
public:
  TAIntegral(){}
  virtual ~TAIntegral(){}

  /// computes integration of f over grid x with n points
  static T Simpson(int n, const double *x, const T *f);
  static T Trapezoid(const FUNCTOR &fun, double a, double b);
  static T Simpson(const FUNCTOR &fun, double a, double b);
  static T Romberg(const FUNCTOR &fun, double a, double b);

protected:
  /// n-stage quadrature of fun over [a,b]. step: h=(b-a)/2^(n-2).
  /// n starts from 0. This method is defined protected, for it's a underlying
  /// method, only supposed to be used incrementally as
  ///  for(i=0;i<m;i++) s=trapzd(fun, a, b, i)
  /// \param FUNTOR should be a class with methods as double operator()(double)
  static T trapzd(const FUNCTOR &fun, double a, double b, int n);
};

#include "TAIntegral.hpp" // definitions of member template functions

#endif
