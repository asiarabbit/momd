/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TATwoPointODE.cxx
  \class TATwoPointODE
  \brief A collection of methods for solving two point value problems. This is a
  abstract base class. Users should supply the ODE and the boundary conditions by
  implementing member function load, score and derivs in the subclass.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/14
  \date Last modified: 2020/09/12 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cmath>
#include <algorithm>
#include "TATwoPointODE.h"
#include "TAException.h"

using std::max;
using std::min;

/// routine for use of TAEqSetSolver::Newton to solve a two point boundary value
/// problem for nvar coupled ODEs by shooting from x1 to x2. Initial values for
/// the nvar ODEs at x2 are generated from the n2 input coefficients v[0..n2-1],
/// using the user-supplied routine score to evaluate the n2 functions f[0..n2-1]
/// that ought to be zero to satisfy the boundary conditions at x2. The functions
/// f are returned on output. TAEqSetSolver::Newton uses a globally convergent
/// Newton's method to adjust the values of v until the function f are zeroed.
/// i.e. TAEqSetSolver::Newton zeros function f(v), which is calculated by routine
/// this routine, namely, shoot. The user-supplied routine derives(x,y,dydx) supplies
/// derivatives information to the ODE integrator
void TATwoPointODE::Shoot(const double *v, double *f){
  fmax = 0; // 0: not store intermediate results
  load(fx1,v,fy); // assign y
  ODEIntegrator(fx1,fx2);
  score(fx2,fy,f); // calculate f
} // end of member function Shoot

void TATwoPointODE::SetRange(double x1, double x2){
  if(x1 == x2) TAException::Error("TATwoPointODE", "SetRange: x1 is equal to x2.");
  fx1 = x1;
  fx2 = x2;
  fh1 = (fx2-fx1)/100.;
} // end of member function SetRange
