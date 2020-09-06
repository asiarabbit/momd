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
/// \param FUNTOR should be a class with methods as double operator()(double)
template<typename T, typename FUNCTOR>
T TAIntegral<T, FUNCTOR>::trapzd(const FUNCTOR &fun, double a, double b, int n){
  if(b - a < 0) TAException::Error("TAIntegral", "trapzd(): b - a < 0");
  if(n < 0) TAException::Error("TAIntegral", "trapzd(): n < 0");

  static T tn; // to store the ultimate result
  if(0 == n) return tn = 0.5*(b-a)*(fun(a)+fun(b)); // the crudest trapezoidal rule

  // T_2n = 1/2*(T_n+H_n), H_n = h*\sum_{i=0}^{n-1}{fun(x_{i+1/2})}
  int nn = 1; // the number of additional interior points: 2 to the power n-2
  for(int i = 0; i < n-1; i++) nn <<= 1; // nn = pow(2, n-2)
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
  static const int npol = 5, nmax = 20;
  static const T epsilon = 1.e-6; // fractional accuracy

  double h[nmax]; T s[nmax]; // h is actually the conventional h^2
  T romb = 0., dromb; // stores the integration result and the error estimate
  h[0] = 1.;
  for(int i = 0; i < nmax; i++){
    s[i] = trapzd(fun, a, b, i);
    if(i >= npol - 1){
      romb = TAInterpolate<T>::PolyInter(h+i-npol, s+i-npol, npol, 0., &dromb);
      if(fabs(dromb) < epsilon*fabs(romb)) return romb;
    }
    // this is the key step. since h^2, not h, is the argument for the extrapolation
    // each time h haves, it entails a factor of (1/2)^2, i.e. 0.25
    h[i+1] = 0.25*h[i];
  } // end while
  TAException::Error("TAIntegral", "Romberg(): Too many steps occurred.");
  return 0.; // never gets here
} // end of the member function Romberg
