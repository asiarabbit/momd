/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAEqSetSolver.hpp
  \class TAEqSetSolver
  \brief Algorithms for root finding and solving nonlinear equation sets. transcribbed
  from Numerical Recipes in C. (Cambridge U. Press, 2002)
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/11
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <algorithm>
#include <cstring>
#include "TAMath.h"
#include "TAException.h"

using std::max;

template <typename FUNC>
TAEqSetSolver<FUNC>::TAEqSetSolver(){}
template <typename FUNC>
TAEqSetSolver<FUNC>::~TAEqSetSolver(){}

/// Given an n-dimensional point x0[0...n-1], the value of the function f0 and
/// the gradient fp0[0...n-1] and the Newton direction p[0...n-1], this method
/// finds a new point x[0...n-1] along the direction p from x0, so that at x
/// the function has decreased "sufficiently". The new fucntion value is returned
/// in f. stepmax is an input quantity that limits the length of the steps so that
/// you do not try to evaluate the function in regions where it is undefined or
/// subject to overflow. The output quantity check is false(0) on a normal exit.
/// It is true(1) when x is too close to x0. In a minimization algorithm,
/// this usually signals convergence and can be ignored. However, in a
/// zero-finding algorithm the calling program should check whether the convergence
/// is spurious. Here f=fmin(x)=1/2*(F.F). where F(x)=0 is the equation set to be solved.
/// NOTE that functor type FUNC should have member function:
/// double operator()(int n, const double *x, double *fvec)
template <typename FUNC>
void TAEqSetSolver<FUNC>::LineSearch(int n, const double *x0, double f0, const double *fp0,
      double *p, double *x, double &f, double stepmax, int &check, FUNC &func){
  // df < alpha* \Nebla f . p = alpha*slope to ensure sufficient decrease in f
  static const double ALF = 1.e-4;
  static const double TOLX = 1.e-7; // convergence criterion on \Delta x.

  check = 0; // initialized by false(0)
  double pnorm = TAMath::norm(n, p);
  if(pnorm > stepmax) for(int i = n; i--;) p[i] *= stepmax / pnorm; // if |p| too big, scale to stepmax
  double slope = TAMath::innerProduct(n, fp0, p); // g'(lambda) = \Nebla f . p
  if(slope >= 0.) TAException::Error("TAEqSetSolver",
    "slope is minus: %f, LineSearch: p not a Newton direction or \
Roundoff errorr problem occurred.", slope);

  double lambda = 0., lambda2, f2, lambda_min, tmpLambda;
  const int nn = n; double fvec[nn]{}; // to store func(x)
  /// calculate lambda_min: lambda_min*p <= x0*TOLX ~1e-7*x0
  double poxMax = 0.; // pox: p over x, p/x
  for(int i = n; i--;) poxMax = max(poxMax, fabs(p[i])/max(fabs(x0[i]), 1.));
  lambda_min = TOLX / poxMax;
  lambda = 1.; // always try Newton step first
  while(1){
    for(int i = n; i--;) x[i] = x0[i] + lambda * p[i]; // x=x0+lambda*p
    f = fmin(n, x, fvec, func);
    if(lambda < lambda_min){ // too small a step, convergence on x reached
      for(int i = n; i--;) x[i] = x0[i];
      check = 1;
      return; // for zero-finding, the calling program should varify the convergence
    } // end if
    else if(f <= f0+ALF*lambda*slope) return; // sufficient function decrease detected
    else{ // find lambda that minimizes f(x0+lambda*p), cached in tmpLambda
      if(1. == lambda) tmpLambda = -slope/(2.*(f-f0-slope)); // the first time for trying the lambda
      else{
        double rhs1 = f -f0-lambda *slope;
        double rhs2 = f2-f0-lambda2*slope;
        double a = (rhs1/(lambda*lambda)-rhs2/(lambda2*lambda2))/(lambda-lambda2);
        double b = (-lambda2*rhs1/(lambda*lambda)+lambda*rhs2/(lambda2*lambda2))/(lambda-lambda2);
        if(0. == a) tmpLambda = -slope/(2.*b);
        else{
          double disc = b*b-3.*a*a*slope;
          if(disc < 0.) tmpLambda = 0.5*lambda; // no extremum in cubic, so take the maximum lambda
          else if(b < 0.) tmpLambda = (-b+sqrt(disc))/(3.*a);
          else tmpLambda = -slope/(b+sqrt(disc)); // numerically more preferrable
        }
        if(tmpLambda > 0.5*lambda) tmpLambda = 0.5*lambda; // lambda <= 0.5*lambda1
      } // end else (a != 0)
    } // end else (lambda finding)
    lambda2 = lambda; f2 = f; // store lambda and f
    lambda = max(tmpLambda, 0.1*lambda); // update lamda with tmpLambda, lambda >= 0.1*lambda1
  } // end while
} // end of the member function LineSearch

/// Given an initial guess x[0...n-1] for a root in n dimensions, find the root
/// of function array func by a globally convergent Newton's method. The output
/// quantity check is false (0) on a normal return and true (1) if the routine
/// has converged to a locla minimum of the function f=func.func. In which case
/// try restarting from a different initial guess.
/// NOTE that functor type FUNC should have member function:
/// double operator()(int n, const double *x, double *fvec)
template <typename FUNC>
void TAEqSetSolver<FUNC>::Newton(double *x, int n, int &check, FUNC &func){
  static const double MAXITS = 200; // maximum iteration times
  static const double TOLF = 1.e-4; // criterion for convergence on func
  static const double TOLMIN = 1.e-6; // check for spurious convergence, i.e. to a local minimum
  static const double TOLX = 1.e-7; // convergence criterion for delta x: |delta x|<TOFLX -> return
  static const double STEPMAX = 100.; // maximum step = STEPMAX * ||x||
  // test for initial guess being a root. Use more stringent test than simply TOLF //
  const int nn = n; double fvec[nn];
  double f = fmin(n, x, fvec, func); // f=1/2*fvec.fvec
  if(TAMath::infNorm(n, fvec) < 0.01*TOLF){ check = 0; return; }

  // start the iterations of the Newton iteration method //
  double stepmax = STEPMAX * max(TAMath::norm(n, x), double(n)); // STEPMAX*|x|
  matrix jj(n,n); // jj(i,j) = partial_func_i/partial_x_j
  double fp[n]{}, x0[n]{}, f0, p[n]{}; // \Nebla fmin, old x and old f; p is for Newton direction
  for(int i = MAXITS; i--;){
    Jacobian(n, x, fvec, jj, func); // n calls of func
    for(int j = 0; j < n; j++){
      for(int k = 0; k < n; k++) fp[j] += jj[k][j] * fvec[k]; // compute \Nebla f
      x0[j] = x[j]; // store x
      p[j] = -fvec[j]; // rhs for Newton Eq: J.dx=-F
    } // end for over j
    f0 = f;
    /// solve the Newton direction p ///
    TAMath::LUSolve(jj, n, p); // Solve J.dx=-F by inversing jj, the solution is stored inplace in p
    LineSearch(n, x0, f0, fp, p, x, f, stepmax, check, func); // f & fvec updated to at new x
    /// test for convergence ///
    if(TAMath::infNorm(n, fvec) < TOLF){ check = 0; return; } // test convergence on fvec
    if(check){ // check for \Nebla_f. If it's ~ zero, possibly a minimum is reached, or it's fine
      double dfmax = 0.; // ||\Nebla_f_i*x_i||_\infty
      for(int j = n; j--;) dfmax = max(dfmax, fabs(fp[j])*max(fabs(x[j]), 1.));
      dfmax /= max(f, 0.5*n); // scale with 1/f, i.e.,
      check = dfmax < TOLMIN ? 1 : 0; // df<TOLMIN: a minimum f; or: possibly a zero
      return;
    } // end the check
    // test convergence on delta_x //
    double dxmax = 0.;
    for(int j = n; j--;) dxmax = max(dxmax, fabs(x[j]-x0[j])/max(fabs(x[j]), 1.));
    if(dxmax < TOLX) return;
  } // end of the loop over iterations
  TAException::Error("TAEqSetSolver", "Newton: Too many iterations, exceeding MAXITS.");
} // end of the member function Newton

/// \retval: 1/2*fvec.fvec
/// \param n: the dimension of the equation set (length of array x and fvec)
/// NOTE that functor type FUNC should have member function:
/// double operator()(int n, const double *x, double *fvec)
template <typename FUNC>
double TAEqSetSolver<FUNC>::fmin(int n, const double *x, double *fvec, FUNC &func){
  func(n, x, fvec); // assign fvec
  return 0.5*TAMath::innerProduct(n, fvec, fvec);
} // end of member function fmin

/// calculate the Jacobian j of function func at x
/// f0 is func at x, and is required to be given upon calling this routine
/// NOTE that functor type FUNC should have member function:
/// double operator()(int n, const double *x, double *fvec)
template <typename FUNC>
void TAEqSetSolver<FUNC>::Jacobian(int n, const double *x, const double *f0, matrix &jj, FUNC &func){
  static const double EPS = 1.e-8; // approximate square root of the machine precision (double)

  const int nn = n;
  double f[nn]{}, xx[nn]{}, h;
  memcpy(xx, x, sizeof(double)*n); // so to keep the const qualifier of double *x
  for(int j = 0; j < n; j++){
    volatile double tmp = xx[j]; // tmp is volatile (readily changeable), so don't optimize tmp
    h = EPS*fabs(tmp);
    if(0. == h) h = EPS;
    xx[j] = tmp + h; // trick to reduce finite precision error
    h = xx[j] - tmp; // so that h is exactly representable in machine
    func(n,xx,f); // calculate f at new x
    xx[j] = tmp;
    for(int i = 0; i < n; i++) jj[i][j] = (f[i] - f0[i]) / h; // forward different formula
  } // end of loop over j
} // end of member function Jacobian
