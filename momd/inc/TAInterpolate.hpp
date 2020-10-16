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
template <typename T>
T TAInterpolate<T>::PolyInter(const double *x, const T *y, int n, double xx, T *dy){
  // find the element in array x that is closest to xx
  int nm = 0;
  double dx, dxm = fabs(xx-x[0]);
  T c[n], d[n]; // c: Pi..(i+m)-Pi(i+m-1), d: Pi..(i+m)-P(i+1)..(i+m)
  for(int i = n; i--;){
    if((dx = fabs(xx-x[i])) < dxm){
      dxm = dx; nm = i;
    } // end for and if
    c[i] = y[i]; d[i] = y[i];
  } // end for over i
  ///// then the Neville's algorithm /////
  // nm-- so as to accommodate Pi..(i+m)=P(i+1)..(i+m)+dm,i
  T result = y[nm--], ddy; // the initial approximation to y; ddy: the correction
  // for m=0, c[i+1]-d[i] reduces to y[i+1]-y[i]
  // so it is ok to initalize c and d to y. c and d starts from m-1
  for(int m = 1; m < n; m++){ // loop over level m
    for(int i = 0; i < n-m; i++){
      double dxi = x[i] - xx, dxmi = x[i+m] - xx;
      T den = dxi - dxmi; // x[i]-x[i+m]
      // two input x's are (to within roundoff) identical
      if(0. == fabs(den)) TAException::Error("TAInterpolate",
        "PolyInter: Elements in array x too close to each other.");
      den = (c[i+1] - d[i]) / den;
      // update c and d to level m
      c[i] = dxi * den;
      d[i] = dxmi * den;
    } // end for over i
    result += ddy = 2*(nm+1) < n-m ? c[nm+1] : d[nm--];
  } // end for over m
  if(dy) *dy = ddy; // the error estimator

  return result;
} // end of member function PolyInter

/// \param len is the length of array x (or y), so that the program would choose
/// the closest interval to envelope xx in the center
/// \param n has the same meaning as in the other overload of PolyInter
template <typename T>
T TAInterpolate<T>::PolyInter(const double *x, const T *y, int len, int n, double xx, T *dy){
  int i = 0, nh = n/2;
  while(x[i++] < xx);
  i = i < nh ? 0 : i - nh;
  if(i + n > len) i = len - n; // so that it does not step out of the border
  return TAInterpolate<T>::PolyInter(x+i, y+i, n, xx, dy);
} // end of member function PolyInter
