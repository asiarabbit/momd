/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAIntegral.hpp
  \class TAIntegral
  \brief Definitions of template member methods of class TAIntegral
  Ref. Numerical Recipes in C 1988-1992, Cambridge U. Press
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/22
  \date Last modified: 2020/07/22 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cmath>
#include "TAInterpolate.h"

/// n-stage quadrature of fun over [a,b]. step: h=(b-a)/2^(n-2).
/// n starts from 0. This method is defined protected, for it's a underlying
/// method, only supposed to be used incrementally as
///  for(i=0;i<m;i++) s=trapzd(fun, a, b, i)
/// \param FUNTOR should be a class with a method: double operator()(double)
template<typename T, typename FUNCTOR>
T TAIntegral<T, FUNCTOR>::trapzd(const FUNCTOR &fun, double a, double b, int n){
  if(b - a < 0) TAException::Error("TAIntegral", "trapzd(): b - a < 0");
  if(n < 0) TAException::Error("TAIntegral", "trapzd(): n < 0");

  static T tn; // to store the ultimate result
  if(0 == n) return tn = 0.5*(b-a)*(fun(a)+fun(b)); // the crudest trapezoidal rule

  // T_2n = 1/2*(T_n+H_n), H_n = h*\sum_{i=0}^{n-1}{fun(x_{i+1/2})}
  int nn = 1; // the number of additional interior points: 2 to the power n-1
  for(int i = n-1; i--;) nn <<= 1; // nn = pow(2, n-1)
  T h = (b-a)/nn, x = a + 0.5*h, hn = 0.;
  for(int i = 0; i < nn; i++, x += h) hn += fun(x);
  tn = 0.5*(tn+h*hn); // Tn->T2n, this replaces s by its refined value
  return tn;
} // end of member function trapzd

// returns the integral of the function fun from a to b with n-stage refinement
// Integraton is performed by the trapezoidal rule
template <typename T, typename FUNCTOR>
T TAIntegral<T, FUNCTOR>::Trapezoid(const FUNCTOR &fun, double a, double b){
  static const double epsilon = 1e-5; // (s_(n+1)-s_(n))/sn: the factional accuracy
  static const int nmax = 20; // maximum number of stages for trapezoidal quadrature

  T tn, tnm = 0.; // Tn and T(n-1)
  int i = 0;
  while(i < nmax){
    tn = trapzd(fun, a, b, i++);
    if(i > 5) // avoid spurious early convergence
      if(fabs(tn-tnm) <= epsilon*tnm) return tn;
    tnm = tn;
  } // end while
  TAException::Error("TAIntegral", "Trapezoid(): Too many steps occurred.");
  return 0.; // never gets here
} // end of member function Trapezoid

// returns the integral of the function fun from a to b with n-stage refinement
// Integraton is performed by the Simpson's rule
template <typename T, typename FUNCTOR>
T TAIntegral<T, FUNCTOR>::Simpson(const FUNCTOR &fun, double a, double b){
  static const T epsilon = 1e-6; // (s_(n+1)-s_(n))/sn: the factional accuracy
  static const int nmax = 20; // maximum number of stages for trapezoidal quadrature

  T tn, sn, tnm = 0., snm = 0.; // Tn and T(n-1), Sn and S(n-1)
  int i = 0;
  while(i < nmax){
    tn = trapzd(fun, a, b, i++);
    sn = (4.*tn-tnm)/3.;
    if(i > 5) // avoid spurious early convergence
      if(fabs(sn-snm) <= epsilon*snm) return sn;
    tnm = tn; snm = sn;
  } // end while
  TAException::Error("TAIntegral", "Simpson(): Too many steps occurred.");
  return 0.; // never gets here
} // end of member function Simpson

template <typename T, typename FUNCTOR>
T TAIntegral<T, FUNCTOR>::Romberg(const FUNCTOR &fun, double a, double b){
  // npol points used in extrapolation and nmax stages in trapezoidal quadrature
  // algebraic precision: 2*np-1, 9 for np = 5
  static const int np = 5, nmax = 20; // np: number of interpolation points
  static const T epsilon = 1.e-6; // fractional accuracy

  double h[nmax]; T s[nmax]; // h is actually the conventional h^2
  T romb = 0., dromb; // stores the integration result and the error estimate
  h[0] = 1.;
  for(int i = 0; i < nmax; i++){
    s[i] = trapzd(fun, a, b, i);
    if(i+1 >= np){
      romb = TAInterpolate<T>::PolyInter(h+i+1-np, s+i+1-np, np, 0., &dromb);
      if(fabs(dromb) < epsilon*fabs(romb)) return romb;
    }
    // this is the key step. since h^2, not h, is the argument for the extrapolation
    // each time h haves, it entails a factor of (1/2)^2, i.e. 0.25
    h[i+1] = 0.25*h[i];
  } // end while
  TAException::Error("TAIntegral", "Romberg(): Too many steps occurred.");
  return 0.; // never gets here
} // end of the member function Romberg

// integral of f(x) over domain [x[0],x[n-1]]. n is the length of x and f
// composed with formula (8.3.6) p.212 Computing Method ver.3 by Guicheng Li
template<typename T, typename FUNCTOR>
T TAIntegral<T, FUNCTOR>::Simpson(int n, const double *x, const T *f){
  // number of individual Simpson intervals
  int k = (n-1) / 2; // f[n-1] is dropped in the case where n is even

  T sum = f[0] + f[2*k];
  for(int i = 0; i < k; i++) sum += 4.*f[2*i+1];
  for(int i = 1; i < k; i++) sum += 2.*f[2*i];

  const double h = (x[n-1] - x[0]) / (n-1);
  if(h < 0) TAException::Error("TAIntegral", "Simpson: x[n-1] is less than x[0].");
  sum *=  h / 3.;
  // Simpson's rule is only for odd n
  // otherwise integral of f[n-1]*h should be explicitly included
  if(n % 2 == 0) sum += f[n-1] * h;

  return sum;
} // end of member function Simpson
