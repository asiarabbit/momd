/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMOMDIS_M.cxx
  \class TAMOMDIS_M
  \brief Calculate core momentum distribution with angular momentum component m
  specified. This is a class to assist class TAMOMDIS.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/10/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cstring>
#include "TAMOMDIS_M.h"
#include "TAMOMDIS.h"
#include "TAFun.h"
#include "TAMath.h"
#include "TABound.h"
#include "TABessel.h"
#include "TASMatrix.h"
#include "TAIntegral.h"
#include "TALegendre.h"

#define bsi TABessel::BesselI
#define sqr TAMath::sqr

static const double INFTY = TABound::INFTY;
static const double PI = TAMath::Pi(), TWOPI = 2.*PI, FOURPI = 4.*PI;
static const cdouble I(0., 1.); // the imaginary number unit

TAMOMDIS_M::TAMOMDIS_M(TAMOMDIS *mom) : fMOM(mom){
  TASMatrix *sc = fMOM->GetSc();
  fRL = sc->GetRL();

  const int ng = fNG = sc->GetNGaus(); // number of gaussians used in the sum-gaus fit
  fAlphajR = new double[ng]; memcpy(fAlphajR, sc->GetAlphajR(), sizeof(double)*ng);
  fAlphajI = new double[ng]; memcpy(fAlphajI, sc->GetAlphajI(), sizeof(double)*ng);
  fBetaM2 = new double[ng];
  for(int j = ng; j--;) fBetaM2[j] = sqr(j/fRL); // 1./(RL/j)^2
  // assign array fMomArr
  const double momMax = 500., h = momMax/kNmom; // MeV/c per nucleon
  fMomArr = new double[kNmom];
  for(int i = kNmom; i--;) fMomArr[i] = h*i; // element with i = 0 stores the total c.s.
} // end of the constructor

TAMOMDIS_M::~TAMOMDIS_M(){
  if(fAlphajR){ delete [] fAlphajR; fAlphajR = nullptr; }
  if(fAlphajI){ delete [] fAlphajI; fAlphajI = nullptr; }
  if(fBetaM2) { delete [] fBetaM2;  fBetaM2  = nullptr; }
} // end of the destructor

double TAMOMDIS_M::ParallelStr(int l, int m, double *momStr){ // accumulate the results
  for(int i = kNmom; i--;) ParallelStr(l, m, fMomArr[i], momStr[i]);
  // or return momStr[0]*2*PI, should be the same
  return TAIntegral<double, int>::Simpson(kNmom, fMomArr, momStr); // the 2nd type arg is irrelevant
} // end of member function ParallelStr
/// the same as ParallelStr, but for diffraction dissociation
double TAMOMDIS_M::ParallelDiff(int l, int m, double *momDiff){
  for(int i = kNmom; i--;) ParallelDiff(l, m, fMomArr[i], momDiff[i]);
  // or return momDiff[0]*2*PI, should be the same
  return TAIntegral<double, int>::Simpson(kNmom, fMomArr, momDiff); // the 2nd type arg is irrelevant
} // end of member function ParallelDiff

/// the integral |B_lmp|^2. Ref[PRC.70.034609-12, the appendix]
inline double besselIp(int m, double x){ // sum of BesslI(m-p, x) over p in [-PMAX, PMAX]
  // modified Bessel function of the first kind, \sum_p{ I_{m-p}(x) }
  static const int PMAX = 7;
  double sum = 0.; // bsi(-m,x) = bsi(m,x)
  for(int p = 1; p <= PMAX; p++) sum += bsi(m+p, x);
  return 2.*sum + bsi(m, x);
} // end of inline function besselIp
// IntToRhoB: so that after the integration, only rho and b are left to be integrated
double TAMOMDIS_M::IntToRhoBStr(double rho, double b, double kz, int l, int m){
  // calculate the coefficients first //
  const double b2 = b*b, rho2 = rho*rho, tworhob = 2.*rho*b;
  cdouble cc(0.,0.);
  for(int j = fNG; j--;){
    double bt2 = fBetaM2[j];
    cc += exp(-b2*bt2)*exp(-rho2*bt2)*besselIp(m,tworhob*bt2)*(fAlphajR[j]+I*fAlphajI[j]);
  } // end for over j
  // then calculate the integral //
  static const auto &intz = [=](double z)->cdouble{
    const double r = sqrt(rho2+z*z);
    return GetRl(r)*TALegendre::Legendre(l,m,z/r)*exp(-I*kz*z);
  }; // end of lambda intz
  return std::norm(cc*TAIntegral<cdouble,decltype(intz)>::Romberg(intz, 0., INFTY));
} // end of member function IntToRhoBStr
// calculate stripping (inelastic) knockout cross sections for momentum kz: dsig/dkz
// calculate stripping (inelastic) knockout cross sections for momentum kz: dsig/dkz
void TAMOMDIS_M::ParallelStr(int l, int m, double kz, double &momStr){
  static const auto &intb = [=](double bn)->double{
    static const auto &intrho = [=](double rho){ return IntToRhoBStr(rho,bn,kz,l,m)*rho; };
    return bn*(1.-std::norm(GetSn(bn)))*
      TAIntegral<double,decltype(intrho)>::Romberg(intrho, 0., INFTY);
  }; // end lambda intb
  momStr = TAMath::Factorial(l-m)/TAMath::Factorial(l+m)/FOURPI* // 2pi/(2l+1)*Clm^2
    TAIntegral<double,decltype(intb)>::Romberg(intb, 0., INFTY);
} // end of member function ParallelStr
// calculate the same result as ParallelStr, but with a somewhat more direct method //
// i.e. sc is integrated over phi without using sum of gaussians to approximate sc
void TAMOMDIS_M::ParallelStr1(int l, int m, double kz, double &momStr){
  // calculate <phi0||sc|^2*|sn|^2|phi0> //
  static const auto &intb = [=](double b)->double{
    static const auto &intrho = [=](double rho)->double{
      const double rho2 = rho*rho, rb2 = rho2 + b*b, tworhob = 2.*rho*b;
      static const auto &intphi = [=](double phi){ // integrate sc over phi
        return std::norm(GetSc(sqrt(rb2-tworhob*cos(phi)))); /// ---> THE DIFFERENCE <----- ///
      }; // end lambda intphi
      static const auto &intz = [=](double z)->cdouble{ // integrate Rl*Plm over z
        const double r = sqrt(rho2+z*z);
        return GetRl(r)*TALegendre::Legendre(l,m,z/r)*exp(-I*kz*z);
      }; // end lambda intz
      return std::norm(TAIntegral<cdouble, decltype(intz)>::Romberg(intz, 0., INFTY))*
        TAIntegral<double, decltype(intphi)>::Romberg(intphi, 2., TWOPI) *rho;
    }; // end lambda intrho
    return b*(1.-std::norm(GetSn(b)))* /// ---> THE DIFFERENCE <----- ///
      TAIntegral<double, decltype(intrho)>::Romberg(intrho, 0., INFTY);
  }; // end lambda intb
  momStr = TAMath::Factorial(l-m)/TAMath::Factorial(l+m)/FOURPI* // 1/(2l+1)*Clm^2
    TAIntegral<double, decltype(intb)>::Romberg(intb, 0., INFTY);
} // end of member function ParallelStr
// calculate stripping (inelastic) knockout cross sections for momentum kz: dsig_diff/dkz
void TAMOMDIS_M::ParallelDiff(int l, int m, double kz, double &momDiff){
  // calculate <phi0||sc|^2*|sn|^2|phi0> //
  static const auto &intb0 = [=](double b)->double{
    static const auto &intrho = [=](double rho)->double{
      double rho2 = rho*rho, b2 = b*b, rb2 = rho2 + b2, tworhob = 2.*rho*b;
      static const auto &intphi = [=](double phi){ // integrate sc over phi
        return std::norm(GetSc(sqrt(rb2-tworhob*cos(phi)))); /// ---> THE DIFFERENCE <----- ///
      }; // end lambda intphi
      static const auto &intz = [=](double z)->cdouble{ // integrate Rl*Plm over z
        const double r = sqrt(rho2+z*z);
        return GetRl(r)*TALegendre::Legendre(l,m,z/r)*exp(-I*kz*z);
      }; // end lambda intz
      return std::norm(TAIntegral<cdouble, decltype(intz)>::Romberg(intz, 0., INFTY))*
        TAIntegral<double, decltype(intphi)>::Romberg(intphi, 2., TWOPI) *rho;
    }; // end lambda intrho
    return b*std::norm(GetSn(b))* /// ---> THE DIFFERENCE <----- ///
      TAIntegral<double, decltype(intrho)>::Romberg(intrho, 0., INFTY);
  }; // end lambda intb0
  double tmp0 = TAMath::Factorial(l-m)/TAMath::Factorial(l+m)/FOURPI* // 1/(2l+1)*Clm^2
    TAIntegral<double, decltype(intb0)>::Romberg(intb0, 0., INFTY);

  // calculate <phi0|sc*sn|phi0> //
  static const auto &intb1 = [=](double b)->cdouble{
    static const auto &intrho = [=](double rho)->cdouble{
      double rho2 = rho*rho, rb2 = rho2 + b*b, tworhob = 2.*rho*b;
      static const auto &intphi = [=](double phi){ // integrate sc over phi
        return GetSc(sqrt(rb2-tworhob*cos(phi))); /// ---> THE DIFFERENCE <----- ///
      }; // end lambda intphi
      static const auto &intz = [=](double z)->cdouble{ // integrate Rl*Plm over z
        const double r = sqrt(rho2+z*z);
        return GetRl(r)*TALegendre::Legendre(l,m,z/r)*exp(-I*kz*z);
      }; // end lambda intz
      return std::norm(TAIntegral<cdouble, decltype(intz)>::Romberg(intz, 0., INFTY))*
        TAIntegral<cdouble, decltype(intphi)>::Romberg(intphi, 2., TWOPI) *rho;
    }; // end lambda intrho
    return b*(GetSn(b))* /// ---> THE DIFFERENCE <----- ///
      TAIntegral<cdouble, decltype(intrho)>::Romberg(intrho, 0., INFTY);
  }; // end lambda intb1
  cdouble tmp1 = TAMath::Factorial(l-m)/TAMath::Factorial(l+m)/FOURPI* // 1/(2l+1)*Clm^2
    TAIntegral<cdouble, decltype(intb1)>::Romberg(intb1, 0., INFTY);

  momDiff = tmp0-std::norm(tmp1); // <phi0||sc|^2|sn|^2|phi0> - |<phi0|sc*sn|phi0>|^2
} // end of member function ParallelDiff

/// the radial wavefunction
double TAMOMDIS_M::GetRl(double r){
  return fMOM->GetRl(r);
} // end of member function GetRl
