/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFNN.cxx
  \class TAFNN
  \brief This class represents nucleon-nucleon scattering amplitude, extracted
  from fitting of pp, pn and nn scattering data over a wide range of energies.
  fNN(q)=kNN/4pi*sigmaNN(i+alphaNN)*exp(-betaNN*q^2)
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/09/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <catch2/catch.hpp>
#include "TAFNN.h"
#include "TAMath.h"
#include "TAException.h"
#include "TAInterpolate.h"

#define polint TAInterpolateD::PolyInter

static const double u = TAMath::uMeV(); // atomic mass unit in MeV
static const int npol = 5; // polynomial interpolation order = npol - 1
static const int npolHalf = npol/2;
static const double fourPi = 4.*TAMath::Pi();

TAFNN::ParOption TAFNN::kParOpt = TAFNN::kHoriuchi;

TAFNN::TAFNN(){}
TAFNN::~TAFNN(){}


// P: projectile, T: target
double TAFNN::IsospinAverage(double vpp, double vpn, int zP, int nP, int zT, int nT){
  double aP = zP + nP, aT = zT + nT;
  if(zP < 0 || zT < 0 || nP < 0 || nT < 0 || aP <= 0 || aT <= 0)
    TAException::Error("TAFNN", "IsospinAverage: One of the following occurs:\
\nzP < 0 || zT < 0 || nP < 0 || nT < 0 || aP <= 0 || aT <= 0");

  return (vpp*(nP*zP+nT*zT) + vpn*(nP*zP+nT*zT)) / (aP*aT);
} // end of member function IsospinAverage

/// isospin-averaged nucleon-nucleon scattering cross section
/// for lab energy pEk.
/// \param isPauli 0: free nn scattering; 1: Pauli-blocking corrected
/// \retval in fm^2
double TAFNN::GetSigmaNN(double pEk, int zP, int nP, int zT, int nT,
    bool isPauli){
  double sigmaPP = 0., sigmaPN = 0.; // pp and np scattering c.s.
  const double b = TAMath::EkPerUToBeta(pEk); // beta=v/c
  const double b2 = b*b, b4 = b2*b2;

  if(!isPauli){
    if(pEk > 10.){
      sigmaPN = -70.67 - 18.18/b + 25.26/b2 + 113.85*b;
      sigmaPP =  13.73 - 15.04/b +  8.76/b2 +  68.67*b4;
    } // end if(pEk > 10.)
    else{
      sigmaPN = 2.73/(0.35*pEk+pow(1-0.0553*pEk,2)) +
        17630./(6.8*pEk+pow(1+0.344*pEk,2));
      sigmaPP =  13.73 - 15.04/b +  8.76/b2 +  68.67*b4;
    } // end else
  } // end if(isPauli)
  else{ // take in Pauli-blocking
    static const double re = 0.17; // effective rho for Pauli-blocking in fm
    sigmaPP =  13.73 - 15.04/b +  8.76/b2 +  68.67*b4;
    sigmaPP *= (1.+7.772*pow(pEk, 0.06)*pow(re, 1.48)) /
      (1.+18.01*pow(re, 1.46));
    sigmaPN = -70.67 - 18.18/b + 25.26/b2 + 113.85*b;
    sigmaPN *= (1.+20.88*pow(pEk, 0.04)*pow(re, 2.02)) /
      (1.+35.86*pow(re, 1.9));
  } /// end the final else
  return IsospinAverage(sigmaPP, sigmaPN, zP, nP, zT, nT) * 0.1;
} // end of member function GetSigmaNN

// various parameterizations exist for alphaNN and betaNN
/// fNN(q)=kNN/4pi*sigmaNN(i+alphaNN)*exp(-betaNN*q^2)
/// alphaNN at arbitrary Ek is interpolated from array alphaNN
double TAFNN::GetAlphaNN(double pEk, int zP, int nP, int zT, int nT){
  if(kHoriuchi == kParOpt) return AlphaBetaHoriuchi(pEk, true);
  return AlphaBetaLenziRay(pEk, zP, nP, zT, nT, true);
} // end of member function GetAlphaNN
double TAFNN::GetBetaNN(double pEk, int zP, int nP, int zT, int nT){
  if(kHoriuchi == kParOpt) return AlphaBetaHoriuchi(pEk, false);
  return AlphaBetaLenziRay(pEk, zP, nP, zT, nT, false);
} // end of member function GetAlphaNN

/// Horiuchi, et al., PRC.75.044607. E: [30, 1000], N=16
double TAFNN::AlphaBetaHoriuchi(double pEk, bool isAlpha){
  if(pEk < 30. || pEk > 1000.) TAException::Error("TAFNN",
    "GetAlphaNN: Ek'd be within [10, 1000] MeV/u for Horiuchi Parameters.");

  static const int n = 16;
  static const double e[n] = { // Elab, in MeV/u
    30., 38., 40., 49., 85., 94., 100., 120., 150., 200.,
    325., 425., 550., 650., 800., 1000.
  };
  static const double a[n] = { // alphaNN, already isospin averaged
    0.87, 0.89, 0.9, 0.94, 1.37, 1.409, 1.435, 1.359, 1.245, 0.953,
    0.305, 0.36, 0.04, -0.095, -0.07, -0.275
  };
  static const double b[n] = { // betaNN, already isospin averaged
    0.685, 0.521, 0.486, 0.390, 0.349, 0.327, 0.322, 0.255, 0.195, 0.125,
    0.075, 0.078, 0.125, 0.16, 0.21, 0.21
  };

  // select a segment to wrap pEk in
  int i = 0; while(e[i++] < pEk); i = i < npolHalf ? 0 : i - npolHalf;
  if(i + npol > n) i = n - npol; // so that it does not step out of the border
  if(isAlpha) return polint(e+i, a+i, npol, pEk);
  return polint(e+i, b+i, npol, pEk);
} // end of member function AlphaBetaHoriuchi

/// Lenzi, et al., PRC.40.2114. E: [10, 94], N=7
/// Ray, et al., PRC.20.1857. E: [100, 2200], N=10
double TAFNN::AlphaBetaLenziRay(double pEk, int zP, int nP, int zT, int nT,
    bool isAlpha){
  if(pEk < 10. || pEk > 2200.) TAException::Error("TAFNN",
    "GetAlphaNN: Ek'd be within [10, 2200] MeV/u for Lenzi & Ray Parameters.");

  static const int n = 17; // 7 + 10
  static const double e[n] = {
    10., 30., 38., 40., 49., 85., 94.,
    100., 150., 200., 325., 425., 550., 650., 800., 1000., 2200.
  };
  static const double app[n] = {
    0.8, 0.87, 0.89, 0.9, 0.94, 1.0, 1.07,
    1.87, 1.53, 1.15, 0.45, 0.47, 0.32, 0.16, 0.06, 0.09, 0.17
  };
  static const double apn[n] = {
    0.8, 0.87, 0.89, 0.9, 0.94, 1.0, 1.07,
    1.00, 0.96, 0.71, 0.16, 0.25, 0.24, 0.35, 0.20, 0.46, 0.50
  };
  static const double bpp[n] = {
    0., 0., 0., 0., 0., 0., 0.,
    0.66, 0.57, 0.56, 0.26, 0.21, 0.04, 0.07, 0.09, 0.09, 0.12
  };
  static const double bpn[n] = {
    0., 0., 0., 0., 0., 0., 0.,
    0.36, 0.58, 0.68, 0.36, 0.27, 0.085, 0.09, 0.12, 0.12, 0.14
  };

  // select a segment to wrap pEk in
  int i = 0; while(e[i++] < pEk); i = i < npolHalf ? 0 : i - npolHalf;
  if(i + npol > n) i = n - npol; // so that it does not step out of the border
  // Lenzi and Ray should be used separately
  if(i <= 6 && i + npol - 1 >= 7){
    if(pEk >= e[7]) i = 7;
    else i = 6 - npol < 0 ? 0 : 6 - npol;
  }
  if(isAlpha) return IsospinAverage(polint(e+i, app+i, npol, pEk),
    polint(e+i, apn+i, npol, pEk), zP, zT, nP, nT);
  return IsospinAverage(polint(e+i, bpp+i, npol, pEk),
    polint(e+i, bpn+i, npol, pEk), zP, zT, nP, nT);
} // end of member function AlphaBetaLenziRay

TEST_CASE("template class TAInterpolate", "[polint]"){
  double x[4] = {1., 2., 3., 4.}, y[4] = {1., 4., 9., 16.}, dy;
  double yy = polint(x, y, 4, -2.5, &dy);
  // TAException::Info("TAInterpolate", "PolyInter: dy is %f", dy);
  CHECK(yy == 6.25);
  CHECK(dy == 0.);
} // end of TEST_CASE
