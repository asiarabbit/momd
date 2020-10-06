/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TABound.h
  \class TABound
  \brief To calculate bound state radial wavefunction of a valence nucleon in
  central potential (Vc+VN+VLS) of a nucleus.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/09/20 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TABound_h_
#define _TABound_h_

#include <string>
#include <vector>
#include "TATwoPointODE.h"

using std::string;
using std::vector;

class TABound : public TATwoPointODE{
public:
  /// inFile: take user input for potentials and single-particle state
  /// of the valence nucleon
  TABound(const vector<double> &vv);
  virtual ~TABound();

  double GetPotential(double r); ///< \retval V(r) = V0+VSO+VC+VL (VL is the centrifugal barrier)
  using TATwoPointODE::operator(); // so that (*this)(int, const double *, double *) is valid
  virtual double operator()(double r){ return GetPotential(r); } ///< interface to minimization routine
  int Getl() const{ return fl; }
  /// \retval Note that Rl=u*r; u = Rl/r. This method returns Rl, i.e. u*r, NOT u
  double GetRl(double r); ///< implement an interpolation

  /// the ODE set: du/dr=u1; du1/dr=-(lambda+v)*u; du2/dr=0,  where E=hbar^2/(2mu)*lambda ///
  /// the two point border condition: u(0)=0; du(0)/dr=1; u(\infty)=0; ///
  /// the ODE set itself: dyi/dx = f_i(x,y1,..yN)
  virtual void derivs(double x, const double *y, double *dydx) override;
  /// calculates the n-vector y[0..n-1] (satisfying the starting boundary conditions, of course)
  /// given the freely specifiable variables of v[0..n2-1] at the initial point x1
  virtual void load(double x1, const double *v, double *y) override;
  /// calculates the discrepancy vector f[0..n2] of the ending boundary conditions,
  /// given the vector y[0..n-1] at the endpoint x2
  virtual void score(double x2, const double *y, double *f) override;

  // solve the radial wavefunction //
  void Bound();

  static const double INFTY; // fx2
private:
  double fZc, fAc; ///< identity of the core
  double fZv, fAv; ///< identity of the valence nucleon
  /// single-particle state for the valence nucleon
  int fn, fl, f2j;
  /// the nuclear potential: WS form
  double fV0, fR0, fA0; ///< central NN force
  double fVS, fRS, fAS; ///< central NN force
  double fRC; ///< charge radius for Coulomb potential

  double fCoeVS, fCoeCoul, fCoel; ///<  VS*RPICOMP^2/AS*<s.l>, Zc*Zv*e2, hbar^2/(2*mu)*l(l+1)
  double f2MuOverHbarSqr; // 2mu/hbar^2

  double fEnergy; ///< the solved eigen energy in MeV
};

#endif
