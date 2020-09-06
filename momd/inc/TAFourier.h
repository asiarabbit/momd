/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFourier.h
  \class TAFourier
  \brief This class implements Fourier transfrom of f(r) to f(q), where f(r) is
  spherically symmetric, and q is in x-y plane.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/21 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAFourier_h_
#define _TAFourier_h_

class TAFourier{
public:
  TAFourier(){}
  virtual ~TAFourier(){}

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
  static double Fourier2D(double q, int nr, const double *r, const double *fr);
  /// put the Fourier transform in an array (fq)
  static void Fourier2D(int nr, const double *r, const double *fr,
    int nq, const double *q, double *fq);
};

#endif
