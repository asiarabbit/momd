/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASpheroidalHarmonics.cxx
  \class TASpheroidalHarmonics
  \brief Usage of two point boundary problems for solving ODE of spheroidal harmonics.
  Currently the shooting method is used.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/14
  \date Last modified: 2020/09/12 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <catch2/catch.hpp>
#include "TASpheroidalHarmonics.h"
#include "TAException.h"
#include "TAEqSetSolver.h"
#include "TAFun.h" // for unit test of TAEqSetSolver::Newton method

using std::cout;
using std::cin;
using std::endl;

// crack the differential equation for spheroidal harmonics using shooting method
// The high-order differential equation are equivalently represented by a coupled
// set of first-order differential equations. Satisfying the boundary condition at
// x1, the rest of the initial conditions are "guessed", then the first-order
//  differential equations are integrated by RK method to x2, so as to satisfy
// the boundary condition at x2. Newton's method are employed to refine the "guess".
// integrate from -1 to 1
void TASpheroidalHarmonics::Sphoot(){
  static const int N2 = 1; // only 1 boundary condition at x2
  int check; // status indicator of Newton's method: 1: possibly a local minimum, not a zero
  double v[N2];
  fdx = 1.e-4;
  SetSave(0, 0.2, 3); // nmax, dxsave, nvar
  SetEPS(1.e-8); // relative error at each step of RKstepper

  fm = 2; fn = 5; fc2 = 16.;
  if(fn < fm || fm < 0){
    TAException::Warn("TASpheroidalHarmonics", "Sphoot: n >= m and m >= 0 are required.");
  }
  fgamma = 1.;
  double q = fn; // compute gamma
  for(int i = 1; i <= fm; i++) fgamma *= -0.5*(fn+i)*q--/i;
  SetRange(-1.+fdx, 0.); // set the inital range: x1, x2
  // this one here is important. To choose an eigenvalue guess properly is an art
  v[0] = fn*(fn+1) - fm*(fm+1) + fc2/2.; // initial guess for eigenvalue: mu = lambda-m(m+1)
  // find v that zeros funtion f in score
  TAEqSetSolver<TATwoPointODE>::Newton(v, N2, check, *this);
  if(check){
    cout << "shoot failed, local minimum instead of zero of f reached." << endl;
    cout << "Bad initial guess for the eigenvalue." << endl;
  }
  fLambda = v[0] + fm*(fm+1.);
} // end of member function Sphoot

/// calculates the discrepancy vector f[0..n2] of the ending boundary conditions,
/// given the vector y[0..n-1] at the endpoint x2
void TASpheroidalHarmonics::score(double x2, const double *y, double *f){
  f[0] = (fn-fm) & 1 ? y[0] : y[1]; // y=0 for odd n-m and y'=0 for even n-m
} // end of member function score
/// calculates the n-vector y[0..n-1] (satisfying the starting boundary conditions, of course)
/// given the freely specifiable variables of v[0..n2-1] at the initial point x1
void TASpheroidalHarmonics::load(double x1, const double *v, double *y){
  // x1 = -1.; v[0] is the eigenvalue, and the guess
  double y0 = (fn-fm) & 1 ? -fgamma : fgamma;
  y[2] = v[0];
  y[1] = -(y[2]-fc2)*y[0]/(2.*(fm+1.));
  y[0] = y0 + y[1]*fdx;
} // end of member function load
/// the ODE set itself: dyi/dx = f_i(x,y1,..yN)
void TASpheroidalHarmonics::derivs(double x, const double *y, double *dydx){
  dydx[0] = y[1];
  dydx[1] = (2.*x*(fm+1.)*y[1]-(y[2]-fc2*x*x)*y[0])/(1.-x*x); // the ode itself
  dydx[2] = 0.;
} // end of member function derivs


TEST_CASE("Newton method for root-finding", "[newt]"){
  class FF : public TAFun<double>{
  public:
    virtual void operator()(int n, const double *x, double *f){
      f[0] = x[0]-exp(x[1]);
      f[1] = x[0]*x[0]+x[1];
    }
  };
  const int n = 2; // dimension of the equation set
  double x[n] = {0., 0.};
  double y[n] = {0., 0.};
  int check;
  FF f;
  TAEqSetSolver<FF>::Newton(x, n, check, f);
  f(n,x,y);
  CHECK(x[0] == Approx(0.652918640419).epsilon(1e-10));
  CHECK(x[1] == Approx(-0.42630275100).epsilon(1e-9));
  CHECK(y[0] == Approx(0.).margin(1.e-7));
  CHECK(y[1] == Approx(0.).margin(1.e-7));
  CHECK(check == 0);
} // end of TEST_CASE[newt]

TEST_CASE("Solve spheroidal harmonics using shoot in TATwoPointODE", "[sphoot]"){
  TASpheroidalHarmonics sp;
  sp.Sphoot();
  CHECK(sp.GetLambda() == Approx(36.9963).margin(0.066));
} // end of TEST_CASE[sphoot]
