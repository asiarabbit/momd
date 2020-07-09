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
  \date Last modified: 2020/07/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TABound.h"

TABound::TABound(const string &inFile){}
virtual ~TABound::TABound(){}

/// construct total potential
void TABound::ConstructPotential(){}
/// calculate and output the bound state radial wavefunction in text file
void TABound::Bound(double *Rl){}
