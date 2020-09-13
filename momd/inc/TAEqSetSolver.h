/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAEqSetSolver.h
  \class TAEqSetSolver
  \brief Algorithms for root finding and solving nonlinear equation sets. Transcribbed
  from Numerical Recipes in C. (Cambridge U. Press, 2002)
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/11
  \date Last modified: 2020/09/03 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAEqSetSolver_h_
#define _TAEqSetSolver_h_

#include "TAMatrix.h"
#include "TATwoPointODE.h"

template <typename FUNC>
class TAEqSetSolver{
public:
  TAEqSetSolver();
  virtual ~TAEqSetSolver();

  /// Given an initial guess x[0...n-1] for a root in n dimensions, find the root
  /// of function array func by a globally convergent Newton's method. The output
  /// quantity check is false (0) on a normal return and true (1) if the routine
  /// has converged to a locla minimum of the function f=func.func. In which case
  /// try restarting from a different initial guess.
  /// NOTE that functor type FUNC should have member function:
  /// double operator()(int n, const double *x, double *kvec)
  static void Newton(double *x, int n, int &check, FUNC &func);
  /// calculate the Jacobian j of function func at x
  /// f0 is func at x, and is required to be given upon calling this routine
  /// NOTE that functor type FUNC should have member function:
  /// double operator()(int n, const double *x, double *kvec)
  static void Jacobian(int n, const double *x, const double *f0, matrix &jj, FUNC &func);

protected:
    /// Given an n-dimensional point x0[0...n-1], the value of the function f0 and
    /// the gradient fp0[0...n-1] and the Newton direction p[0...n-1], this method
    /// finds a new point x[0...n-1] along the direction p from x0, so that at x
    /// the function has decreased "sufficiently". The new fucntion value is returned
    /// in f. stepmax is an input quantity that limits the length of the steps so that
    /// you do not try to evaluate the function in regions where it is undefined or
    /// subject to overflow. The output quantity check is false (0) on a normal exit.
    /// It is true (1) when x is too close to x0. In a minimization algorithm,
    /// this usually signals convergence and can be ignored. However, in  a
    /// zero-finding algorithm the calling program should check whether the convergence
    /// is spurious. Here f=fmin(x)=1/2*(F.F). where F(x)=0 is the equation set to be solved.
    /// NOTE that functor type FUNC should have member function:
    /// double operator()(int n, const double *x, double *kvec)
    static void LineSearch(int n, const double *x0, double f0, const double *fp0,
          double *p, double *x, double &f, double stepmax, int &check, FUNC &func);

  /// \retval: 1/2*kvec.kvec
  /// \param n: the dimension of the equation set (length of array x and kvec)
  /// NOTE that functor type FUNC should have member function:
  /// double operator()(int n, const double *x, double *kvec)
  static double fmin(int n, const double *x, double *kvec, FUNC &func);

protected:
  static double *kvec; // to pass function arrays assigned in func
};

#include "TAEqSetSolver.hpp" // the implementations

#endif
