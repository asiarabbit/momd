/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TATwoPointODE.h
  \class TATwoPointODE
  \brief A collection of methods for solving two point value problems. This is a
  abstract base class. Users should supply the ODE and the boundary conditions by
  implementing member function load, score and derivs in the subclass.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/14
  \date Last modified: 2020/08/14 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TATwoPointODE_h_
#define _TATwoPointODE_h_

#include "TAODESolver.h"

class TATwoPointODE : public TAODESolver{
public:
  TATwoPointODE() : fx1(-9999.), fx2(-9999.){}
  virtual ~TATwoPointODE(){}

  /// routine for use of TAEqSetSolver::Newton to solve a two point boundary value
  /// problem for nvar coupled ODEs by shooting from x1 to x2. Initial values for
  /// the nvar ODEs at x2 are generated from the n2 input coefficients v[0..n2-1],
  /// using the user-supplied routine score to evaluate the n2 functions f[0..n2-1]
  /// that ought to be zero to satisfy the boundary conditions at x2. The functions
  /// f are returned on output. TAEqSetSolver::Newton uses a globally convergent
  /// Newton's method to adjust the values of v until the function f are zeroed.
  /// i.e. TAEqSetSolver::Newton zeros function f(v), which is calculated by routine
  /// this routine, namely, shoot. The user-supplied routine derives(x,y,dydx) supplies
  /// derivatives information to the ODE integrator.
  void Shoot(int n, const double *v, double *f);
  void operator()(int n, const double *v, double *f){ Shoot(n, v, f); }

  /// the following are user-specific routines dependent of the certain ODEs to be solved ///
  /// calculates the discrepancy vector f[0..n2] of the ending boundary conditions,
  /// given the vector y[0..n-1] at the endpoint x2
  virtual void score(double x2, const double *y, double *f) = 0;

protected:
  /// integration range
  double fx1, fx2;
};

#endif
