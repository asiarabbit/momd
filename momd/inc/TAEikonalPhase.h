/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAEikonalPhase.h
  \class TAEikonalPhase
  \brief compute eikonal phase \chi=\chi_N + \chi_C, S=exp(i*\chi), where \chi_N
  represents nuclear phase, and \chi_C Coulomb phase.
  This is a class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/25
  \date Last modified: 2020/09/05 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAEikonalPhase_h_
#define _TAEikonalPhase_h_

#include <complex>

typedef std::complex<double> cdouble;

class TAEikonalPhase{
public:
  TAEikonalPhase(int zP, int aP, int zT, int aT, double ek);
  virtual ~TAEikonalPhase();

  /// \retval eikonal nucclear phase = 1/KNN*\int_0^\infty{q*rhoP*rhoT*fNN*J0(qb)}
  /// \param b: impact parameter
  virtual cdouble GetPhaseN(double b);
  /// \retval eikonal Coulumb phase = 2\eta*ln(kb)
  /// \param b: impact parameter
  virtual double GetPhaseC(double b);
  /// \retval \chi_C+\chi_N
  virtual cdouble GetPhase(double b);
  // Fourier transform of nucleus density on the transverse plane
  // alphap is the width (sigma) of proton (nucleon) Gaussian density
  void SetProjectileDenstiy(int n, const double *r, const double *rho);
  void SetTargetDenstiy(int n, const double *r, const double *rho);
  void SetNucleonSize(double alphap){ fAlphaP = alphap; }

  friend class TAPhaseN;

protected:
  int fZP, fAP, fZT, fAT; // projectile and target (N, Z)
  double fEk; ///< Elab in MeV/u
  double fMu; // the reduced mass in MeV (mP*mT)/(mP+mT)
  /// Array R in those densities are supposed to be equidistant
  int fNRP, fNRT; ///< length of fRP and fRT
  double *fRP, *fRhoP; ///< density of the projectile (fRP, fRhoP) in r (fm)
  double *fRT, *fRhoT; ///< density of the target (fRT, fRhoT) in r (fm)
  double fAlphaP; // the width (sigma) of proton Gaussian distribution

  double fAlphaNN, fBetaNN, fSigmaNN;
  double fBeta; // the relative speed: P w.r.t. T
};

#endif
