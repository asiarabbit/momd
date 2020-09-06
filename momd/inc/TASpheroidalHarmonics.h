/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASpheroidalHarmonics.h
  \class TASpheroidalHarmonics
  \brief Usage of two point boundary problems for solving ODE of speroidal harmonics:
    (1-x^2)*d^2y/dx^2-2(m+1)x*dy/dx+(mu-c^2*x^2)y=0
  Currently the shooting method is used.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/14
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASperoidalHarmonics_h_
#define _TASperoidalHarmonics_h_

#include "TATwoPointODE.h"

class TASperoidalHarmonics : public TATwoPointODE{
public:
  TASperoidalHarmonics(){}
  virtual ~TASperoidalHarmonics(){}

  void Sphoot();
  /// calculates the discrepancy vector f[0..n2] of the ending boundary conditions,
  /// given the vector y[0..n-1] at the endpoint x2
  void score(double x2, const double *y, double *f);
  /// calculates the n-vector y[0..n-1] (satisfying the starting boundary conditions, of course)
  /// given the freely specifiable variables of v[0..n2-1] at the initial point x1
  void load(double x1, const double *v, double *y);
  /// the ODE set itself: dyi/dx = f_i(x,y1,..yN)
  void derivs(double x, const double *y, double *dydx);

protected:
  /// boundary condition at +-1 is actually applied at a distance dx from the boundary
  double fdx;
  /// m and c2 (c^2) are equation parameters. n is the quantum number of the eigenvalue
  int fm, fn;
  double fc2;
  double fgamma; // y(1) are normaolized to gamma = (-)^m*(n+m)!/(2^m*m!*(n-m)!)
};

#endif
