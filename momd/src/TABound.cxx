/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TABound.cxx
  \class TABound
  \brief To calculate bound state radial wavefunction of a valence nucleon in
  central potential (Vc+VN+VLS) of a nucleus.    This is a standalone class. It
  is supposed to solicit potential inputs and valence nucleon state quanta from
  users and output the results in arrays to text files.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TABound.h"

TABound::TABound(const vector<double> &vv) : fV(0), fRl(0){
  // (c+v) is P (projectile) //
  fZc = vv[5]; fAc = vv[6];
  fZv = vv[0] - fZc; fAv = vv[1] - fAc;
  // the quantum state of the valence nucleon //
  fn = vv[8]; fl = vv[9]; f2j = vv[10];
  // assign the potential parameters //
  fV0 = vv[11]; fR0 = vv[12]; fA0 = vv[13];
  fVS = vv[14]; fRS = vv[15]; fAS = vv[16];
  fRC = vv[17];
} // end of the constructor
TABound::~TABound(){
  if(fV){ delete [] fV; fV = nullptr; }
  if(fRl){ delete [] fRl; fRl = nullptr; }
}

/// construct total potential in fV
void TABound::ConstructPotential(){
  if(fV) return;
  fV = new double[kNV];
} // end of member function ConstructPotential
/// calculate and output the bound state radial wavefunction in fRl
void TABound::Bound(){
  if(fRl) return;
  fRl = new double[kNRl];

  ConstructPotential();
} // end of member function
/// output the solved wavefunction
const double *TABound::GetRl(){
  if(!fRl) Bound();
  return fRl;
} // end of member function GetRl
