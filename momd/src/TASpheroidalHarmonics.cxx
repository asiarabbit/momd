/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASpheroidalHarmonics.cxx
  \class TASpheroidalHarmonics
  \brief Usage of two point boundary problems for solving ODE of speroidal harmonics.
  Currently the shooting method is used.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/14
  \date Last modified: 2020/08/14 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include "TASpheroidalHarmonics.h"
#include "TAException.h"
#include "TAEqSetSolver.h"

using std::cout;
using std::cin;
using std::endl;

// crack the differential equation for speroidal harmonics using shooting method
// The high-order differential equation are equivalently represented by a coupled
// set of first-order differential equations. Satisfying the boundary condition at
// x1, the rest of the initial conditions are "guessed", then the first-order
//  differential equations are integrated by RK method to x2, so as to satisfy
// the boundary condition at x2. Newton's method are employed to refine the "guess".
// integrate from -1 to 1
void TASperoidalHarmonics::Sphoot(){
  static const int N2 = 1; // only 1 boundary condition at x2
  int check; // status indicator of Newton's method: 1: possibly a local minimum, not a zero
  double v[N2];
  fdx = 1.e-4; fnvar = 3;

  while(1){
    cout << "Input m, n, and c-squared: " << endl;
    if(!(cin >> fm >> fn >> fc2)) break;
    if(fn < fm || fm < 0) continue;
    fgamma = 1.; double q = fn; // compute gamma
    for(int i = 1; i <= fm; i++) fgamma *= -0.5*(fn+i)*(q--/i);
    fx1 = -1.+fdx; fx2 = 0.; // ste the inital range
    v[0] = fn*(fn+1) - fm*(fm+1) + fc2/2.; // initial guess for eigenvalue: mu = lambda-m(m+1)
    // find v that zeros funtion f in score
    TAEqSetSolver<TATwoPointODE>::Newton(v, N2, check, *this);
    if(check){
      cout << "shoot failed, local minimum instead of zero of f reached." << endl;
      cout << "Bad initial guess for the eigenvalue." << endl;
    }
    else cout << "mu(m,n): " << v[0] << endl;
  } // end while
} // end of member function Sphoot

/// calculates the discrepancy vector f[0..n2] of the ending boundary conditions,
/// given the vector y[0..n-1] at the endpoint x2
void TASperoidalHarmonics::score(double x2, const double *y, double *f){
  f[0] = (fn-fm) & 1 ? y[0] : y[1];
} // end of member function score
/// calculates the n-vector y[0..n-1] (satisfying the starting boundary conditions, of course)
/// given the freely specifiable variables of v[0..n2-1] at the initial point x1
void TASperoidalHarmonics::load(double x1, const double *v, double *y){
  double y0 = (fn-fm) & 1 ? -fgamma : fgamma;
  y[2] = v[0];
  y[1] = (y[2]-fc2)/(2.*(fm+1.))*y[0];
  y[0] = y0 + y[1]*fdx;
} // end of member function load
/// the ODE set itself: dyi/dx = f_i(x,y1,..yN)
void TASperoidalHarmonics::derivs(double x, const double *y, double *dydx){
  dydx[0] = y[1];
  dydx[1] = (2.*x*(fm+1.)*y[1]-(y[2]-fc2*x*x)*y[0])/(1.-x*x);
  dydx[0] = 0.;
} // end of member function derivs
