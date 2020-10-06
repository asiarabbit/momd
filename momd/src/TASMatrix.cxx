/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASMatrix.cxx
  \class TASMatrix
  \brief This class calculate the eikonal S-matrix. It takes nucleus densities
  and nucleon-nucleon scattering amplitude (fNN) as input. The output is stored
  in arrays.   Since the S-matrix is fitted with sum of Gaussians, this class
  takes care of the fitting as well.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/10/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <fstream>
#include <cstring>
#include "TASMatrix.h"
#include "TAEikonalPhase.h"
#include "TAFit.h"
#include "TABound.h"

static const cdouble I(0., 1.); // imaginary number unit I

TASMatrix::TASMatrix(int zP, int aP, int zT, int aT, double ek){
  fNG = 21; fRL = 50.;
  fNB = 200; fBMax = TABound::INFTY; // fm
  fEikonalPhase = new TAEikonalPhase(zP, aP, zT, aT, ek);
} // end of the constructor
TASMatrix::~TASMatrix(){
  if(fEikonalPhase){ delete fEikonalPhase; fEikonalPhase = nullptr; }
  if(fAlphajR){ delete [] fAlphajR; fAlphajR = nullptr; }
  if(fAlphajI){ delete [] fAlphajI; fAlphajI = nullptr; }
} // end of the destructor

/// calculate and output the S-matrix
cdouble TASMatrix::SMatrix(double b) const{
  return exp(I*fEikonalPhase->GetPhase(b));
} // end of member function SMatrix
///< real and imag part of S-matrix, in an array (0, kBMax) with step h=kBMax/fNB
void TASMatrix::SMatrix(double *b, double *smr, double *smi) const{
  const double h = fBMax / fNB;
  double bb = 0.;
  for(int i = 0; i < fNB; i++){
    cdouble sm = SMatrix(b[i]=bb);
    smr[i] = sm.real();
    smi[i] = sm.imag();
    bb += h;
  } // end for over i
} // end of member function SMatrix(double*, double*)
/// fit the S-matrix to an expansion of Gaussians
/// the results are stored in member arrays alphaj (a complex number)
void TASMatrix::GaussianFit(){
  static bool isCalled = false;
  if(isCalled) return;

  // calculate a list of S-Matrix w.r.t. b (impact parameter)
  // real and imag part of the s-matrix
  const int n = fNB, m = fNG; // nb: bins of the impact parameter; ng: number of base functions (gaussians)
  if(!fAlphajR && !fAlphajI){
    fAlphajR = new double[m]; fAlphajI = new double[m];
  }
  else{
    TAException::Warn("TASMatrix", "Gaussian: fAlphaj has been assigned.");
    return;
  }
  bool isFit[m]; memset(isFit, 1, sizeof(isFit)); // whether it is for fit, or fixed
  // alphaj for the 1-st gaussian is fixed to 1., as usually S-Matrix approaches 1. for large b
  isFit[0] = false; fAlphajR[0] = 1.; fAlphajI[0] = 0.;
  matrix covarR(m, m), covarI(m, m);
  double b[n], smr[n], smi[n], dy[n]; // dy: the variance of y (smr or smi)
  SMatrix(b, smr, smi); // calculate the S-Matrix, assign b, smr and smi
  double chi2R = -9999., chi2I = -9999.;
  for(int i = 1; i <= n; i++) dy[i] = 1.;
  TAFit::LSMFit(b, smr, dy, n, m, fAlphajR, isFit, covarR, chi2R, gaus, &fRL);
  TAFit::LSMFit(b, smi, dy, n, m, fAlphajI, isFit, covarI, chi2I, gaus, &fRL);

  // now we get alpha in fAlphajR and fAlphajI, and we are done here
  isCalled = true;
} // end of member function GaussianFit

double *TASMatrix::GetAlphajR(){
  if(!fAlphajR) GaussianFit();
  return fAlphajR;
} // end of member function GetAlphajR
double *TASMatrix::GetAlphajI(){
  if(!fAlphajI) GaussianFit();
  return fAlphajI;
} // end of member function GetAlphajI

/// gaussians used for S-Matrix expansion, \retval array funci[ma] evaluated at x
void TASMatrix::gaus(double x, double *funci, int ma, const double *p){
  for(int i = 0; i < ma; i++) funci[i] = exp(-x*x*i/p[0]);
} // end of member function gaus
