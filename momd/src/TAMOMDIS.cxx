/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMOMDIS.cxx
  \class TAMOMDIS
  \brief This is a global class, to take control of the whole flow of the
  program. It is also responsible for generating various momentum distributions.
  This is supposed to be a singleton class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/07/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TAMOMDIS.h"
#include "TAOutput.h"

TAMOMDIS *TAMOMIS::kInstance = nullptr;

TAMOMDIS *TAMOMDIS::Instance(){
  if(!kInstance) kInstance = new TAMOMDIS();
  return kInstance;
}

virtual ~TAMOMDIS::TAMOMDIS(){}

void TAMOMDIS::Go(){
  MOMParallel();
}

////------------ PREPARE FOR THE MOMDIS INTEGRAL ----------------/////
void TAMOMDIS::Prepare(){
  static bool isCalled = false;
  if(isCalled) return;
  // to extract the radial wavefunction
  TABound *fBound = new Bound(fConfigDir);
  fRl = new double[TABound::kNRl];
  fBound->Bound(fRl);

  // to calculate the S-matrix and fit it for alphaj and betaj
  TASMatrix *fSMatrix = new TASMatrix(fConfigDir);
  fAlphaj = new double[TASMatrix::kNG];
  fBetaj = new double[TASMatrix::kNG];
  GaussianFit(fAlphaj, fBetaj);

  const int fl = fBound->Getl();

  isCalled = true;
} // end of member function Prepare

///--- CALCULATE PARTIAL C.S. && CORE MOMENTUM DISTRIBUTION w.r.t. m ---///
void TAMOMIS::Parallel(){
  Prepare(); // calculate Rl and S-matrix

  // calculate momdis and sigma for each m
  TAMOMDIS_M *mom = new TAMOMDIS_M(this);
  // the array length of momentum distribution
  const int kNmom = TAMOMIS_M::kNmom;
  double sigmaStr_M[fl+1]{}, sigmaDiff_M[fl+1]{}; // +-m -> 2*m
  double momStr_M[fl+1][kNmom]{}, momDiff_M[fl+1][kNmom]{}; // +-m -> 2*m
  double momTotal_M[fl+1][kNmom]{}; // sum of stripping and diffraction momdis
  double sigmaStr = 0., sigmaDiff = 0., sigmaTotal = 0.;
  // accumulated over m
  double momStr[kNmom]{}, momDiff[kNmom]{}, momTotal[kNmom]{};
  for(int m = 0; m <= l; m++){
    sigmaStr_M[m] = mom->ParallelStr(m, momStr_M[m]);
    sigmaDiff_M[m] = mom->ParallelDiff(m, momDiff_M[m]);
    for(int i = kNmom; i--;) momTotal_M[m][i] = momStr_M[m][i] + momDiff_M[m][i];
    if(0 != m){ // take in the -m part
      sigmaStr_M[m] *= 2.;
      sigmaDiff_M[m] *= 2.;
      sigmaTotal_M[m] *= 2.;
      for(int i = kNmom; i--;) momTotal_M[m][i] *= 2.;
    } // end if
    // accumulate over m
    sigmaStr += sigmaStr_M[m];
    sigmaDiff += sigmaDiff_M[m];
    sigmaTotal += sigmaTotal_M[m];
    for(int i = kNmom; i--;){
      momStr[i] += momStr_M[m][i];
      momDiff[i] += momDiff_M[m][i];
      momTotal[i] += momTotal_M[m][i];
    } // end for over i
  } // end for over m

  /////----- OUTPUT THE RESULT -----/////
  TAOutput::PrintKOCS(fl, sigmaStr, sigmaStr_M, sigmaDiff_M,
    sigmaDiff,sigmaTotal_M, sigmaTotal);
  TAOutput::PrintToFile(kNmom, mom->GetMomArr(), momTotal);
} // end of the member function Parallel
