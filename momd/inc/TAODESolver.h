/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAODESolver.h
  \class TAODESolver
  \brief A collection of methods for solving ordinary differential equations
  based on Runge-Kutta methods. Note that this class is dedicated for initial
  value problems. For two point boundary value problems, please resort to class
  TATwoPointODE.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/08
  \date Last modified: 2020/08/14 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAODESolver_h_
#define _TAODESolver_h_

#include "TAMatrix.h"

class TAODESolver{
public:
  TAODESolver();
  virtual ~TAODESolver();

  /// a stepper program for RK method, transcribbed from routine rkqs in p.719 of
  /// Numerical Recipes in C
  /// Given n dependent variables y[0...n-1] and their derivatives dydx[0...n-1] at x,
  /// this method step y by hdid (initial step is htry, step suggested for next
  /// call of rk is hnext), and update x as x+=hdid.
  /// Thi method features 'quality-controlled' adaptive stepsize control.
  /// eps is the error tolerance w.r.t. max_i{yerr[i]/yscale[i]}. yerr is provided
  /// by another routine of rk algorithm.
  /// \param eps is a constant fractional error: eps = err/(y+dy), for RK5, eps~O(h^5)
  void RKStepper(double *y, double *dydx, int n, double &x, double htry,
      double eps, const double *yscale, double &hdid, double &hnext);
  /// Given values for n variables y[0...n-1] and their derivatives dydx[0...n-1]
  /// at x, this method advances the solution over an interval h and return the
  /// incremented variables as yout with the fifth-order Cash-Karp Runge-Kutta method,
  /// and and estimated local truncation error yerr usinig the embedded
  /// fourth-order RK method. The user supplies the routine derivs(x,y,dydx) to
  /// provide dy/dx at x.
  void RKCashKarp(const double *y, double *dydx, int n, double &x, double h, double *yout, double *yerr);

  /// Runge-Kutta driver with adaptive stepsize control. Integrate starting values
  /// ystart[0...n-1] from x1 to x2 with constant fractional accuracy eps, storing intermediate
  /// results in memeber vector kxp and matrix kyp;
  /// \param eps is a constant fractional error: eps = err/(y+dy), for RK5, eps~O(h^5)
  /// \param h1 the initial step, hmin: the minimum step allowed
  /// \param nOK and nBAD: counts of steps where input h is adopted (hdid==h)
  /// this method uses RKStepper for implementation of the RK method
  /// solution at x2 would be updated inplace in ystart
  /// \param kmax = 0 would cancel the storage of the intermediate results
  void ODEIntegrator(double x1, double x2){
    ODEIntegrator(fy,fnvar,x1,x2,EPS,fh1,fhmin,fnOK,fnBAD);
  }
  void ODEIntegrator(double *ystart, int nvar, double x1, double x2, double eps,
      double h1, double hmin, int &nOK, int &nBAD);
  /// \retval output the mapping x->y[ndim], and returns the length of the arrays
  int GetSolution(double *x, double *y, int ndim = 0) const;
  void Solve(double x1, double x2, const double *v); // v: the inital conditions

  /// \brief set parameters involving sainvg the intermediate results of the ODE solution
  /// \param nvar number of coupled ODEs
  /// If kmax != 0 results are stored at approximate intervals kdxsave within kmax steps
  /// illegal values passed would be ignored
  void SetSave(int maxStepCount, double maxStepSize, int nvar);
  void SetEPS(double eps); ///< relative error for each step
  void SetMinStep(double hmin); ///< minimum step allowed
  const double *GetYatX2(){ return fy; }

protected:
  /// the ODE set itself: dyi/dx = f_i(x,y1,..yN)
  virtual void derivs(double x, const double *y, double *dydx) = 0;
  /// calculates the n-vector y[0..n-1] (satisfying the starting boundary conditions, of course)
  /// given the freely specifiable variables of v[0..n2-1] at the initial point x1
  virtual void load(double x1, const double *v, double *y) = 0;
  double *fy; ///< for storing the solution vector
  double EPS; ///< relative error for each step
  double fh1, fhmin; // fhmin: allowd minimum step size
  int fnOK, fnBAD; // nOK and nBAD: counts of steps where input h is adopted and rejected

  // for output the resutls of the solution to the ODEs //
  /// Preset kmax and kdxsave in the calling program. If kmax != 0 results are stored
  /// at approximate intervals kdxsave in the arrays kxp[0..kcount-1] and yp[0..knvar-1][0..kcount],
  /// where kcount is output by member function ODEIntegrator
  ///!!! value of kmax and kdxsave should be specified in the calling program !!!///
  int fmax; ///< maximum length of the solution array to be stored
  double fdxsave; ///< only saves solutions for steps wtih h > dxsave, user specific
  int fnvar; // the dimenstion of the ODE
  int fcount; /// real dimension of the solution array (fcount*fnvar)
  ///!!! memory allocation of fxp and fyp should be in the calling program !!!///
  double *fxp; ///< for storing intermediate results of the solution of the ODEs
  /// a kmax*n matrix, for storing the intermediate results of the dependent variables
  matrix *fyp; ///< to store the ODE integration results, dim: [0..knvar-1][0..kmax]
};

#endif
