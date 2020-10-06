/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMOMDIS.cxx
  \class TAMOMDIS
  \brief This is a global class, to take control of the whole flow of the
  program. It is also responsible for generating various momentum distributions.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/10/05 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>
#include "TAMOMDIS.h"
#include "TAMOMDIS_M.h"
#include "TASMatrix.h"
#include "TABound.h"
#include "TAOutput.h"
#include "TAException.h"

using std::ifstream;
using std::vector;
using std::istringstream;

TAMOMDIS::TAMOMDIS(const string &configFile) : fConfigFile(configFile), fBound(0),
    fSc(0), fSn(0), fMOM_M(0){
  fVecConfig.reserve(19); // at least 19 parameters would be pushed into this vector
} // end of the constructor

TAMOMDIS::~TAMOMDIS(){
  if(fSc){ delete fSc; fSc = nullptr; }
  if(fSn){ delete fSn; fSn = nullptr; }
  if(fBound){ delete fBound; fBound = nullptr; }
  if(fMOM_M){ delete fMOM_M; fMOM_M = nullptr; }
} // end of the destructor

void TAMOMDIS::Go(){
  Parallel();
} // end of member function Go

////------------ PREPARE FOR THE MOMDIS INTEGRAL ----------------/////
// solve the boundwave function for the valence nucleon and calculate the S-Matrix
void TAMOMDIS::Configure(){
  static bool isCalled = false;
  if(isCalled) return;

  LoadConfigFile(); // assign fVecConfig
  // to calculate the S-matrix and fit it for alphaj and betaj //
  // vv[0-4]: ZP, AP, ZT, AT, Ek, Zc, Ac
  const auto &vv = fVecConfig;
  const int zP = vv[0], aP = vv[1], zT = vv[2], aT = vv[3], zc = vv[5], ac = vv[6];
  const double Ek = vv[4]; // MeV/u
  const int zv = zP - zc, av = aP - ac;
  fSc = new TASMatrix(zc, ac, zT, aT, Ek);
  fSn = new TASMatrix(zv, av, zT, aT, Ek);
  // to solve the radial wavefunction
  fBound = new TABound(vv);
  fMOM_M = new TAMOMDIS_M(this);

  isCalled = true;
} // end of member function Prepare

///--- CALCULATE PARTIAL C.S. && CORE MOMENTUM DISTRIBUTION w.r.t. m ---///
void TAMOMDIS::Parallel(){
  Configure(); // calculate Rl and S-matrix
  int l = fBound->Getl(); // the orbit angular momentum quantum number

  // calculate momdis and sigma for each m
  // the array length of momentum distribution
  const int kNmom = TAMOMDIS_M::kNmom;
  // -m and m have the same dsigma/dkz
  double sigmaStr_M[l+1]{}, sigmaDiff_M[l+1]{}, sigmaTotal_M[l+1]{};
  double momStr_M[l+1][kNmom]{}, momDiff_M[l+1][kNmom]{};
  double momTotal_M[l+1][kNmom]{}; // sum of stripping and diffraction momdis
  double sigmaStr = 0., sigmaDiff = 0., sigmaTotal = 0.;
  // accumulate over m
  double momStr[kNmom]{}, momDiff[kNmom]{}, momTotal[kNmom]{};
  for(int m = 0; m <= l; m++){
    sigmaStr_M[m] = fMOM_M->ParallelStr(l, m, momStr_M[m]);
    sigmaDiff_M[m] = fMOM_M->ParallelDiff(l, m, momDiff_M[m]);
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
  TAOutput::PrintKOCS(l, sigmaStr_M, sigmaStr, sigmaDiff_M,
    sigmaDiff, sigmaTotal_M, sigmaTotal);
  TAOutput::PrintToFile(kNmom, fMOM_M->GetMomArr(), momTotal, "mom.txt");
} // end of the member function Parallel

// read configFile to fVecConfig //
// skip spaces and tabs, return subscript of the 1st valid char
inline int skipCrap(const string &s){
  int tmp = 0; char c;
  for(;;c = s[tmp++]) if(' ' != c && '\t' != c) break;
  return tmp - 1;
}
void TAMOMDIS::LoadConfigFile(){
  static bool isCalled = false;
  if(isCalled) return;

  ifstream fin(fConfigFile.c_str());
  string line;
  double item;
  while(getline(fin, line)){
    if(!line.size()) continue;
    const char c = line[skipCrap(line)];
    if('#' == c || '\0' == c) continue;

    istringstream sline(line);
    while(sline >> item) fVecConfig.push_back(item);
  } // end while

  if(!fVecConfig.size())
    TAException::Warn("TAMOMDIS", "LoadConfigFile: fVecConfig size is 0. Empty config file?");
  isCalled = true;
} // end of the memeber function LoadConfigFile

TABound *TAMOMDIS::GetBound(){
  if(!fBound) Configure();
  return fBound;
} // end of member function GetBound
double TAMOMDIS::GetRl(double r){ return GetBound()->GetRl(r); } ///< the radial wavefunction
/// returns the S-matrix of the core-target
TASMatrix *TAMOMDIS::GetSc(){
  if(!fSc) Configure();
  return fSc;
} // end of member function GetSMatrix
/// returns the S-matrix of the valence nucleon-target
TASMatrix *TAMOMDIS::GetSn(){
  if(!fSn) Configure();
  return fSn;
} // end of member function GetSMatrix
cdouble TAMOMDIS::GetSc(double b){ return GetSc()->SMatrix(b); }
cdouble TAMOMDIS::GetSn(double b){ return GetSn()->SMatrix(b); }
