/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TABound.cxx
  \class TABound
  \brief To calculate bound state radial wavefunction of a valence nucleon in
  central potential (Vc+VN+VLS) of a nucleus.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/10/05 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TABound.h"
#include "TAMath.h"
#include "TAMinimize.h"
#include "TAInterpolate.h"
#include "TAEqSetSolver.h"

#define polint TAInterpolateD::PolyInter
#define sqr TAMath::sqr

// physical constants for calculation of the potential
static const double hbarc = TAMath::hbarc();
static const double rp2 = sqr(hbarc/TAMath::PionPMMass()); // pion compton wavelength^2
static const double e2 = hbarc*TAMath::FineStructureConstant(); // e^2/(4pi*epsilon0)

const double TABound::INFTY = 200.; // in fm, the x2, physical infinity to nuclear systems

TABound::TABound(const vector<double> &vv){
  // (c+v) is P (projectile) //
  fZc = vv[5]; fAc = vv[6];
  fZv = vv[0] - fZc; fAv = vv[1] - fAc;
  // the quantum state of the valence nucleon //
  fn = vv[8]; fl = vv[9]; f2j = vv[10];
  // assign the potential parameters //
  fV0 = vv[11]; fR0 = vv[12]; fA0 = vv[13];
  fVS = vv[14]; fRS = vv[15]; fAS = vv[16];
  fRC = vv[17];

  double j = f2j/2.;
  fCoeVS = fVS*rp2/fAS * (j*(j+1.)-fl*(fl+1.)-0.75)/2.;
  fCoeCoul = fZc*fZv*e2; // VS*RPICOMP^2/AS, Zc*Zv*e2;
  const double mu = (fAc*fAv)/(fAc+fAv)*TAMath::uMeV(); // the reduced mass in MeV
  f2MuOverHbarSqr = 2.*mu/sqr(hbarc);
  fCoel = fl*(fl+1.);
} // end of the constructor
TABound::~TABound(){}

///< \retval V(r) = V0+VSO+VC+VL (VL is the centrifugal barrier)
double TABound::GetPotential(double r){
  const double v0 = fV0*TAMath::WoodsSaxon(r, fA0, fR0); // the central NN potential
  // the spin-orbin potential //
  double vs = exp((r-fRS)/fAS); // vs = -vs0*r^{pi}_{compton}^2/r*dws/dr
  vs = -fCoeVS * vs/(sqr(1.+vs)*r);
  // the Coulomb potential //
  double vc = fCoeCoul;
  if(r >= fRC) vc /= r; else vc *= (3.-sqr(r/fRC)) / (2.*fRC);
  // the centrifugal barrier //
  double vl = fCoel/(r*r);

  return (v0 + vs + vc)*f2MuOverHbarSqr + vl;
} // end of member function GetPotential

// ---------- the defintion of the eigenvalue problem ---------------- //
/// the ODE set: du/dr=u1; du1/dr=-(lambda+v)*u; du2/dr=0,  where E=hbar^2/(2mu)*lambda ///
/// the two point border condition: u(0)=0; du(0)/dr=1; u(\infty)=0; ///
/// the ODE set itself: dyi/dx = f_i(x,y1,..yN)
void TABound::derivs(double x, const double *y, double *dydx){
  dydx[0] = y[1];
  dydx[1] = -(y[2]+GetPotential(x)) * y[0];
  dydx[2] = 0.;
} // end of member function derivs
/// calculates the n-vector y[0..n-1] (satisfying the starting boundary conditions, of course)
/// given the freely specifiable variables of v[0..n2-1] at the initial point x1
/// x1 = 0., x2 = \infty = 200.
void TABound::load(double x1, const double *v, double *y){
  y[0] = v[0];
  y[1] = 1.;
  y[2] = v[1];
} // end of member function load
/// calculates the discrepancy vector f[0..n2] of the ending boundary conditions,
/// given the vector y[0..n-1] at the endpoint x2
void TABound::score(double x2, const double *y, double *f){
  f[0] = y[0];
} // end of member function score

void TABound::Bound(){
  static const int N2 = 1; // only 1 boundary condition at x2
  int check; // status indicator of Newton's method: 1: possibly a local minimum, not a zero
  SetNvar(3); // nvar, number of ODEs in the ODE set
  SetEPS(1.e-8); // relative error at each step of RKstepper
  SetRange(0., INFTY);
  double rmin; // where the potential is at its minimum
  // the eigenvalue, lambda
  double v[N2] = {TAMinimize<TABound>::Brent(fx1,fx2,*this,rmin)};

  TAEqSetSolver<TABound>::Newton(v, N2, check, *this);
  if(check)
    TAException::Info("TABound", "Bound: shoot failed, local minimum instead of\
zero of f reached.\n Bad initial guess for the eigenvalue.");

  // output the results //
  SetNmax(1000); // fy has been assigned
  ODEIntegrator(fx1, fx2); // the wavefunction is now available via GetSolution()
  fEnergy = v[0]*f2MuOverHbarSqr; // and the eigen energy
} // end of member function

/// \retval Note that Rl=u*r; u = Rl/r. This method returns Rl, i.e. u*r, NOT u
/// implement an interpolation
double TABound::GetRl(double r){ // u*r
  static const int NPOL = 4; //  // polynomial interpolation order = npol - 1: 4 is enough!
  if(r < 0.) TAException::Error("TABound", "GetRl: input r(%f) is minus.", r);
  if(r > INFTY) TAException::Error("TABound", "GetRl: input r(%f) exceeds INFTY.", r);
  return polint(fxp, (*fyp)[0].at(0), fcount, NPOL, r)*r;
} // end of member function GetRl
