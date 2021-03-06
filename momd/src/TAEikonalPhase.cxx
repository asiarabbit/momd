/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAEikonalPhase.cxx
  \class TAEikonalPhase
  \brief compute eikonal phase \chi=\chi_N + \chi_C, S=exp(i*\chi), where \chi_N
  represents nuclear phase, and \chi_C Coulomb phase.
  This is a static class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/25
  \date Last modified: 2020/09/10 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cstring>
#include <catch2/catch.hpp>
#include "TAFun.h"
#include "TAEikonalPhase.h"
#include "TAFNN.h"
#include "TABessel.h"
#include "TAIntegral.h"
#include "TAMath.h"
#include "TAFourier.h"

// fine structure constant ~1/137, e^2/hbarc
static const double alpha = TAMath::FineStructureConstant();
static const double hbarc = TAMath::hbarc();
static const double u = TAMath::uMeV();
static const double fourPi = 4.*TAMath::Pi();
static const cdouble I(0., 1.); // imaginary number unit I

TAEikonalPhase::TAEikonalPhase(int zP, int aP, int zT, int aT, double ek)
    : fZP(zP), fAP(aP), fZT(zT), fAT(aT), fEk(ek){
  fMu = fAP*fAT/(fAP+fAT)*u;
  if(zP < 0 || aP < 0 || zT < 0 || aT < 0 || ek < 0. || aP <= zP || aT <= zT)
    TAException::Error("TAEikonalPhase",
      "TAEikonalPhase: abnormaol input.");

  int nP = aP - zP, nT = aT - zT;
  // fNN(q)=kNN/4pi*sigmaNN(i+alphaNN)*exp(-betaNN*q^2)
  fAlphaNN = TAFNN::GetAlphaNN(ek, zP, nP, zT, nT);
  fBetaNN  =  TAFNN::GetBetaNN(ek, zP, nP, zT, nT);
  fSigmaNN = TAFNN::GetSigmaNN(ek, zP, nP, zT, nT);
  fBeta = TAMath::EkPerUToBeta(ek);

  fNRP = fNRT = 0;
  fRP = fRhoP = fRT = fRhoT = nullptr;
} // end of the constructor

TAEikonalPhase::~TAEikonalPhase(){
  if(fRP){ delete [] fRP; fRP = nullptr; }
  if(fRhoP){ delete [] fRhoP; fRhoP = nullptr; }
  if(fRT){ delete [] fRT; fRT = nullptr; }
  if(fRhoT){ delete [] fRhoT; fRhoT = nullptr; }
} // end of the destructor

void TAEikonalPhase::SetProjectileDenstiy(int n, const double *r, const double *rho){
  if(n <= 0) TAException::Error("TAEikonalPhase",
    "SetProjectileDenstiy: array length n <= 0, n = %d", n);
  if(fRP){
    TAException::Warn("TAEikonalPhase", "SetProjectileDenstiy: fRP is not nullptr");
    delete [] fRT; fRT = nullptr;
  }
  if(fRhoP){
    TAException::Warn("TAEikonalPhase", "SetProjectileDenstiy: fRhoP is not nullptr");
    delete [] fRhoP; fRhoP = nullptr;
  }
  fNRP = n;
  fRP = new double[n]{};
  fRhoP = new double[n]{};
  memcpy(fRP, r, n*sizeof(double));
  memcpy(fRhoP, rho, n*sizeof(double));
  // normalize to nucleon number
  const double cc = fAP/TAFourier::Fourier2D(0, fNRP, fRP, fRhoP);
  for(int i = n; i--;) fRhoP[i] *= cc;
} // end of member function SetProjectileDenstiy
void TAEikonalPhase::SetTargetDenstiy(int n, const double *r, const double *rho){
  if(n <= 0) TAException::Error("TAEikonalPhase",
    "SetTargetDenstiy: array length n <= 0, n = %d", n);
  if(fRT){
    TAException::Warn("TAEikonalPhase", "SetTargetDenstiy: fRT is not nullptr");
    delete [] fRT; fRT = nullptr;
  }
  if(fRhoT){
    TAException::Warn("TAEikonalPhase", "SetTargetDenstiy: fRhoT is not nullptr");
    delete [] fRhoT; fRhoT = nullptr;
  }
  fNRT = n;
  fRT = new double[n]{};
  fRhoT = new double[n]{};
  memcpy(fRT, r, n*sizeof(double));
  memcpy(fRhoT, rho, n*sizeof(double));
  // normalize to nucleon number
  const double cc = fAT/TAFourier::Fourier2D(0, fNRT, fRT, fRhoT);
  for(int i = n; i--;) fRhoT[i] *= cc;
} // end of member function SetTargetDenstiy

/// \retval eikonal nucclear phase = 1/KNN*\int_0^\infty{q*rhoP*rhoT*fNN*J0(qb)}
/// \param b: impact parameter
cdouble TAEikonalPhase::GetPhaseN(double b){
  if(!fNRP || !fNRT) TAException::Error("TAEikonalPhase",
    "GetPhaseN: Densities may not be assigned yet.");
  // Romberg: compute integral \int_0^\infty{q*rhoP*rhoT*fNN*J0(qb)}
  static const auto &ph = [=](double q)->double{
    const double pden = exp(-TAMath::sqr(q*fAlphaP)/4.); // nucleon size factor
    return q*
      TAFourier::Fourier2D(q, fNRP, fRP, fRhoP)*pden * // rhoPq
      TAFourier::Fourier2D(q, fNRT, fRT, fRhoT)*pden * // rhoTq
      exp(-fBetaNN*q*q)*TABessel::BesselJ0(q*b);
  }; // end of the lambda expression
  return TAIntegral<double, decltype(ph)>::Romberg(ph, 0., 300.) *
    fSigmaNN/fourPi*(I+fAlphaNN);
} // end member function GetPhaseN

/// \retval eikonal Coulumb phase = 2\eta*ln(kb)
/// \param b: impact parameter
double TAEikonalPhase::GetPhaseC(double b){
  // Sommerfeld parameter eta = zp*zt*e^2/(hbar*v)=zp*zt/beta * alpha
  double eta = fZP*fZT/fBeta * alpha; // alpha: fine structure constant
  double k = fMu*fBeta/hbarc; // wave number of AP-AT scattering system, in fm^(-1)
  return 2.*eta*log(k*b); //2.*\eta*ln(k*b)
} // end of member function GetPhaseC

/// \retval \chi_C+\chi_N
/// \param b: impact parameter
cdouble TAEikonalPhase:: GetPhase(double b){
  return GetPhaseN(b) + GetPhaseC(b);
} // end of memeber function GetPhase

TEST_CASE("template class TAIntegral", "[quad]"){
  SECTION("polynomial tests"){
    class QuadTest : public TAFun<double>{
    public:
      double operator()(double x) const override{
        static const int n = 9;
        return (n+1)*pow(x, n);
      } // end of member function
    }; // end of the class
    QuadTest quad;
    double t = TAIntegral<double, QuadTest>::Trapezoid(quad, 0., 1.);
    double s = TAIntegral<double, QuadTest>::Simpson(quad, 0., 1.);
    double r = TAIntegral<double, QuadTest>::Romberg(quad, 0., 1.);
    CHECK(t == Approx(1.).epsilon(1.e-5));
    CHECK(s == Approx(1.).epsilon(1.e-6));
    CHECK(r == 1.);
  } // end section0
  SECTION("trigonometric function tests"){
    class QuadTest : public TAFun<double>{
    public:
      double operator()(double x) const override{
        return sin(x);
      } // end of member function
    }; // end of the class
    QuadTest quad;
    const double b = TAMath::Pi()/2.;
    double t = TAIntegral<double, QuadTest>::Trapezoid(quad, 0., b);
    double s = TAIntegral<double, QuadTest>::Simpson(quad, 0., b);
    double r = TAIntegral<double, QuadTest>::Romberg(quad, 0., b);
    CHECK(t == Approx(1.).epsilon(1.e-5));
    CHECK(s == Approx(1.).epsilon(1.e-6));
    CHECK(r == Approx(1.).epsilon(1.e-6));
  } // end section0
} // end of TEST_CASE
