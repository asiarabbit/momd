/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAInterpolate.h
  \class TAInterpolate
  \brief Polynomial interpolation. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/07/19 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAInterpolate_h_
#define _TAInterpolate_h_

#include <complex>

template <typename T>
class TAInterpolate{
public:
  TAInterpolate(){}
  virtual ~TAInterpolate(){}

  /// Polynomial interpolation using Neville's algorithm
  /// given n points in arrays x and y, and input x, this routine returns the
  /// interpolated func value y at xx, and assigning the error estimate to dy.
  /// Ref. Numerical Receipes in C: p109
  static T PolyInter(const double *x, const T *y, int n, double xx, T *dy = nullptr);
  /// \param len is the length of array x (or y), so that the program would choose
  /// the closest interval to envelope xx in the center
  static T PolyInter(const double *x, const T *y, int len, int n, double xx, T *dy = nullptr);
};

#include "TAInterpolate.hpp"

typedef TAInterpolate<double> TAInterpolateD;
typedef TAInterpolate<std::complex<double>> TAInterpolateC;

#endif
