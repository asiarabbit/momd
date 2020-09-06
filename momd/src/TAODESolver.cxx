/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAODESolver.cxx
  \class TAODESolver
  \brief A collection of methods for solving ordinary differential equations
  based on Runge-Kutta methods.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/08
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cmath>
#include <algorithm>
#include "TAODESolver.h"
#include "TAException.h"

using std::max;
using std::min;

TAODESolver::TAODESolver() : fxp(nullptr), fyp(nullptr){
  fmax = 1000; fdxsave = 0.2; fnvar = 0;
} // end of the constructor
TAODESolver::~TAODESolver(){
  if(fxp){ delete [] fxp; fxp = nullptr; }
  if(fyp){ delete    fyp; fyp = nullptr; }
} // end of the destrcutor


/// a stepper program for RK method, transcribbed from routine rkqs in p.719 of
/// Numerical Recipes in C
/// Given n dependent variables y[0...n-1] and their derivatives dydx[0...n-1] at x,
/// this method step y by hdid (initial step is htry, step suggested for next
/// call of rk is hnext), and update x as x+=hdid.
/// Thi method features 'quality-controlled' adaptive stepsize control.
/// eps is the error tolerance w.r.t. max_i{yerr[i]/yscale[i]}. yerr is provided
/// by another routine of rk algorithm.
/// \param eps is a fractional error: eps = err/(y+dy), for RK5, eps~O(h^5)
void TAODESolver::RKStepper(double *y, double *dydx, int n, double &x, double htry,
    double eps, double *yscale, double &hdid, double &hnext){
  static const double SAFETY = 0.9; // to scale h so as to be smaller, and safer
  static const double EXPL = -0.2, EXPH = -0.25;
  static const double ERRLIM = 1.889e-4; // =pow(5/SAFETY, 1/EXPH), to control h not too inflated

  double errmax, h = htry; // set stepsize to the initial trial value
  const int nn = n; double ytemp[nn], yerr[nn]; // new y and its truncation error

  // complete one step with adaptive stepsize //
  while(1){
    RKCashKarp(y, dydx, n, x, h, ytemp, yerr); // take a step
    errmax = 0.;
    // find the "worst-offender" in array yerr
    for(int i = 1; i < n; i++) errmax = max(errmax, fabs(yerr[i]/yscale[i]));
    errmax /= eps; // Delta1/Delta0 = yerr/(y+dy) / epsilon  for constant fractional error
    if(errmax <= 1.) break;
    double htemp = SAFETY*h*pow(errmax, EXPH);
    h = h >= 0. ? max(htemp, 0.1*h) : min(htemp, 0.1*h); // no more than a factor of 10.
    if(double(x + h) == x) TAException::Error("TAODESolver", "RKStepper: Stepsize underflow.");
  } // end while
  // update h -- so long as while is exit, h is to be relaxed (get bigger)
  if(errmax > ERRLIM) hnext = SAFETY*h*pow(errmax, EXPL);
  else hnext = 5.*h; // no more than a factor of 5
  x += hdid = h; // update x
  for(int i = n; i--;) y[i] = ytemp[i]; // update y
} // end of member function RKStepper

/// Given values for n variables y[0...n-1] and their derivatives dydx[0...n-1]
/// at x, this method advances the solution over an interval h and return the
/// incremented variables as yout with the fifth-order Cash-Karp Runge-Kutta method,
/// and and estimated local truncation error yerr usinig the embedded
/// fourth-order RK method. The user supplies the routine derivs(x,y,dydx) to
/// provide dy/dx at x.
/// \param eps is a constant fractional error: eps = err/(y+dy), for RK5, eps~O(h^5)
void TAODESolver::RKCashKarp(double *y, double *dydx, int n, double &x, double h,
    double *yout, double *yerr){
  // Cash-Karp Parameters for Embedded Runge-Kutta Method //
  // k1=h*f(xn, yn); k2=h*f(xn+a2*h, yn+b21*k1) ... k6=h*f(xn+a6*h,yn+b61*k1+...+b65*k5)
  // y_(n+1)=yn+c1*k1+c2*k2+...+c6*k6, y'_(n+1)=yn+c1'*k1+c2'*k2+...+c6'*k6,
  static const double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1., a6 = 0.875;
  static const double b21 = 0.2;
  static const double b31 = 3./40., b32 = 9./40.;
  static const double b41 = 0.3, b42 = -0.9, b43 = 1.2;
  static const double b51 = -11./54., b52 = 2.5, b53 = -70./27., b54 = 35./27.;
  static const double b61 = -1631./55296., b62 = 175./512., b63 = 575./13824.,
    b64 = 44275./110592., b65 = 253./4096.;
  static const double c1 = 37./378., c3 = 250./621, c4 = 125./594., c6 = 512./1771.;
  // the error estimate: \Delta = y_(n+1)-y'_(n+1) = \sum_i{(c_i-c'_i)*k_i}=\sum_i{dc_i*k_i}
  static const double dc1 = c1-2825./27648., dc3 = c3-18575./48384.,
    dc4 = c4-13525./55296., dc5 = -277./14336., dc6 = c6-0.25;

  const int nn = n;
  // note that k2-6 here are just derivatives, without h, different from the formula
  double k2[nn]{}, k3[nn]{}, k4[nn]{}, k5[nn]{}, k6[nn]{}, ytemp[nn]{};
  // five steps to accomplish fifth-order Cash-Karp RK method //
  // the first step
  for(int i = 0; i < n; i++) ytemp[i] = y[i]+h* b21*dydx[i];
  derivs(x+a2*h, ytemp, k2);
  // the second step
  for(int i = 0; i < n; i++) ytemp[i] = y[i]+h*(b31*dydx[i]+b32*k2[i]);
  derivs(x+a3*h, ytemp, k3);
  // the third step
  for(int i = 0; i < n; i++) ytemp[i] = y[i]+h*(b41*dydx[i]+b42*k2[i]+b43*k3[i]);
  derivs(x+a4*h, ytemp, k4);
  // the fourth step
  for(int i = 0; i < n; i++) ytemp[i] = y[i]+h*(b51*dydx[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]);
  derivs(x+a5*h, ytemp, k5);
  // the fifth step
  for(int i = 0; i < n; i++) ytemp[i] = y[i]+h*(b61*dydx[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i]);
  derivs(x+a6*h, ytemp, k6);

  for(int i = 0; i < n; i++){
    // accumulates increments with proper weights
    yout[i] = y[i] + h*(c1*dydx[i]+c3*k3[i]+c4*k4[i]+c6*k6[i]);
    // estimate error as difference between fourth and fifth order methods
    yerr[i] = h*(dc1*dydx[i]+dc3*k3[i]+dc4*k4[i]+dc5*k5[i]+dc6*k6[i]);
  } // end for over i
} // end of member function RKCashKarp

/// Runge-Kutta driver with adaptive stepsize control. Integrate starting values
/// ystart[0...n-1] from x1 to x2 with constant fractional accuracy eps, storing intermediate
/// results in static memeber vector fxp and matrix fyp;
/// \param eps is a constant fractional error: eps = err/(y+dy), for RK5, eps~O(h^5)
/// \param h1 the initial step, hmin: the minimum step allowed
/// \param nOK and nBAD: counts of steps where input h is adopted (hdid==h)
/// this method uses stepper for implementation of the RK (or other) method and
/// the corresponding stepsize control solution at x2 would be updated inplace in ystart
void TAODESolver::ODEIntegrator(double *ystart, int nvar, double x1, double x2, double eps,
    double h1, double hmin, int &nOK, int &nBAD){
  static const int MAXSTEP = 10000, TINY = 1.e-30; // maximum steps allowed over the solution

  const int n = nvar; // count of the dependent variables
  SetSave(fmax, fdxsave, nvar); // allocate memory to fxp and fyp
  double y[n], yscale[n], dydx[n]; // dy/dx, yscale is to scale y[i] for a uniform error estimate
  for(int i = 0; i < n; i++) y[i] = ystart[i];
  double x = x1, h = x2 - x1 > 0. ? fabs(h1) : -fabs(h1), hdid, hnext;
  // for saving the results //
  fyp = new matrix(fmax, n);
  fxp = new double[fmax];
  double xsave; // x that can be saved
  if(fmax > 0) xsave = x1 - fdxsave*2.; // so that x1 can be saved for sure
  nOK = nBAD = fcount = 0;
  // here we starts the stepping loop //
  for(int i = 0; i < MAXSTEP; i++){ // take at most MAXSTEP steps
    derivs(x, y, dydx); // called in the first place to obtain yscale
    for(int j = 0; j < n; j++) yscale[j] = fabs(y[j]) + fabs(h*dydx[j]) + TINY; // FIXME
    // store the intermediate results //
    if(fmax > 0 && fcount < fmax && fabs(x-xsave) > fabs(fdxsave)){
      fxp[fcount] = x;
      for(int j = 0; j < n; j++) (*fyp)[j][fcount] = y[i];
      fcount++;
      xsave = x;
    } // end if
    if((x+h-x2)*(x+h-x1) > 0.) h = x2-x; // if stepsize can overshoot, decrease
    RKStepper(y, dydx, n, x, h, eps, yscale, hdid, hnext);
    if(hdid == h) nOK++; else nBAD++;
    if((x-x2)*(x2-x1) >= 0.){ // then the end of the interval has been reached
      for(int j = 0; j < n; j++) ystart[i] = y[i];
      if(fmax > 0){ // save the final step before ending the program
        fxp[fcount] = x;
        for(int j = 0; j < n; j++) (*fyp)[j][fcount] = y[i];
        fcount++;
      } // end if
      return;
    } // end if
    if(fabs(hnext) <= hmin)
      TAException::Error("TAODESolver", "ODEIntegrator: stepsize too small.");
    h = hnext;
  } // end for over i
  TAException::Error("TAODESolver", "ODEIntegrator: Too many steps occurred.");
} // end of member function ODEIntegrator

int TAODESolver::GetSolution(double *x, double *y, int ndim) const{
  if(fcount <= 0) TAException::Error("TAODESolver",
    "GetSolution: Length of the solution array is zero or minus: %d", fcount);
  if(ndim < 0 || ndim >= fnvar) TAException::Error("TAODESolver",
    "GetSolution: abnormal ndim, should be within [0, fnvar), ndim: %d", ndim);
  for(int i = fcount; i--;){
    x[i] = fxp[i];
    y[i] = (*fyp)[ndim][i];
  } // end for over i
  return fcount;
} // end of member function GetSolution

void TAODESolver::SetSave(int maxStepCount, int maxStepSize, int nvar){
  if(maxStepSize != fdxsave) fdxsave = maxStepSize;
  if(fnvar == nvar && maxStepCount == fmax) return; // no need for memoroy change

  fnvar = nvar; fmax = maxStepCount;
  if(fxp){ delete [] fxp; fxp = new double[fmax]{}; }
  if(fyp){ delete    fyp; fyp = new matrix(nvar, fmax); }
} // end of member function SetSave


/// Initial values for the nvar ODEs at x2 are generated from the n2 input coefficients
/// v[0..n2-1] by load(...). The user-supplied routine derives(x,y,dydx) supplies
/// derivatives information to the ODE integrator
void TAODESolver::Solve(double x1, double x2, const double *v){
  static const double EPS = 1.e-6;

  if(0 == fnvar) TAException::Error("TAODESolver", "Solve: fnvar not assigned.");
  if(x1 == x2) TAException::Error("TAODESolver", "Solve: x1 is equal to x2.");
  int nbad, nok; // nOK and nBAD: counts of steps where input h is adopted
  double h1 = (x2-x1)/100., hmin = 0., y[fnvar];

  fmax = 0; // 0: not store intermediate results
  load(x1,v,y); // assign y
  ODEIntegrator(y,fnvar,x1,x2,EPS,h1,hmin,nok,nbad);
} // end of member function Shoot
