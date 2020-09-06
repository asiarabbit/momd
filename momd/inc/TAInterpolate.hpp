/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAInterpolate.hpp
  \class TAInterpolate
  \brief Polynomial interpolation. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/09/03 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cmath>
#include "TAMath.h"
#include "TAException.h"

/// Polynomial interpolation using Neville's algorithm
/// given n points in arrays x and y, and input x, this routine returns the
/// interpolated func value y at xx, and assigning the error estimate to dy.
/// Ref. Numerical Receipes in C: p109
template<typename T>
T TAInterpolate<T>::PolyInter(const double *x, const T *y, int n,
    double xx, T *dy){
  // find the element in array x that is closest to xx
  int nm = 0; double dx, dxm = fabs(xx-x[0]);
  T c[n], d[n]; // c: Pi..(i+m)-Pi(i+m-1), d: Pi..(i+m)-P(i+1)..(i+m)
  for(int i = n; i--;){
    if((dx = fabs(x[i]-xx)) < dxm){
      dxm = dx; nm = i;
    } // end for and if
    c[i] = y[i]; d[i] = y[i];
  } // end for over i
  ///// then the Neville's algorithm /////
  // nm-- so as to accommodate Pi..(i+m)=P(i+1)..(i+m)+dm,i
  T result = y[nm--], ddy; // the initial approximation to y; ddy: the correction
  // for m=0, c[i+1]-d[i] reduces to y[i+1]-y[i]
  // so it is ok to initalize c and d to y
  for(int m = 0; m < n-1; m++){ // loop over level m
    for(int i = 0; i < n-m-1; i++){ // so that (i+m)_max = n-1
      double dxi = x[i] - xx, dxmi = x[i+m] - xx;
      T den = dxi - dxmi; // x[i]-x[i+m]
      // two input x's are (to within roundoff) identical
      if(0. == fabs(den)) TAException::Error("TAInterplate",
        "PolyInter: Elements in array x too close to each other.");
      den = (c[i+1] - d[i]) / den;
      // update c and d to level m
      c[i] = dxi * den;
      d[i] = dxmi * den;
    } // end for over i
    ddy = 2*nm < n-m ? c[nm+1] : d[nm--];
    result += ddy;
  } // end for over m
  if(dy) *dy = ddy; // the error estimator

  return result;
} // end of member function PolyInter
