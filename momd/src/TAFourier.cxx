/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFourier.cxx
  \class TAFourier
  \brief This class implements Fourier transfrom of f(r) to f(q), where f(r) is
  spherically symmetric, and q is in x-y plane.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/25 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TAFourier.h"
#include "TAMath.h"
#include "TAException.h"

static const double fourPi = 4.*TAMath::Pi();

// computes the Fourier transform of a spherically symmetric function fr,
// defined on a grid of r with nr points. Returns the Fourier  transform in
// function fq on a grid q with nq points.
// the transform is of the form:
// f(q) = \int[dr\vec f(r) exp(-iq\vec.s\vec)] = 4pi/q\int[dr rf(r)sin(qr)]
// integration is over 0->\infty
//
// implements Filon's Integration formula in M Abramowitz and I A Stegun (1970)
// Handbook of Mathematical Functions, Dover, New York, p894
// accuracy~O(h^3) for odd nr (i.e., 0->2n), which drops to O(h^2) for even nr
// where the last term is added using rectangle quadrature.
// so we suggest that even nr be avoided.
double TAFourier::Fourier2D(double q, int nr, const double *r, const double *fr){
  const double h = r[1] - r[0]; // r and k step
  if(h < 0.)
    TAException::Error("TAFourier", "Fourier2D: dr is minus.\
The r array must be of ascending order.");
  for(int i = 1; i < nr-1; i++) if(r[i+1] - r[i] != h)
    TAException::Error("TAFourier", "Fourier2D: array r is not equidistant.");

  int n = (nr-1)/2; // nr: 0->2n (+1)
  const double k = q; // just a different notation
  double kh = k*h, kh2 = kh*kh, kh3 = kh2*kh, fq = 0.; // fq: stores the final result
  double sinkh = sin(kh), coskh = cos(kh);
  double alpha = 1./kh + sinkh*coskh/kh2 - 2.*sinkh*sinkh/kh3;
  double beta  = (1+coskh*coskh)/kh2 - 2.*sinkh*coskh/kh3; // beta/2 in Filon
  double gamma = 4.*(sinkh/kh3 - coskh/kh2);
  double sEven = 0, sOdd = 0.; // f(r)sin(kr): 0,2,...,2n; 1,3,..2n-1
  int je, jo; // for even and odd times of h
  for(int j = 1; j <= n; j++){ // j==0 omitted where sin(kr)=0
    je = 2*j; jo = je - 1;
    sEven += je*h*fr[je]*sin(je*kh); // rf(r)sinkr
    sOdd  += jo*h*fr[jo]*sin(jo*kh); // rf(r)sinkr
  } // end for over r
  double rfrR = je*h*fr[je];
  sEven = -0.5*rfrR*sin(je*kh);
  fq = fourPi/k* h*(-alpha*rfrR*cos(kh*2.*n) + beta*sEven + gamma*sOdd);
  // complete the integral with the last interval: Rf(R)sin(kR)*h
  if(nr%2 == 0) fq += fourPi/k* h*(nr-1)*fr[nr-1]*sin(kh*(nr-1)) *h;
  return fq;
} // end of member function Fourier2D(...)

/// put the Fourier transform in an array
void TAFourier::Fourier2D(int nr, const double *r, const double *fr, int nq,
    const double *q, double *fq){
  const double h = r[3] - r[2], dk = q[3] - q[2]; // r and k step
  if(h < 0. || dk < 0.)
    TAException::Error("TAFourier", "Fourier2D(...): dr or dq is minus.\
The r and q array must be of ascending order.");

  // calculate the integral for each momentum k
  for(int i = 0; i < nq; i++) fq[i] = Fourier2D(q[i], nr, r, fr);
} // end of member function Fourier2D(...)
