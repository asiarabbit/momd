/**
  SUNNY Project, Anyang Normal University, IMP-CAS
	\file TAManyBodySD.C
	\class TAManyBodySD
	\brief Slater determinant (SD) class for many-body problems. Each SD represents
	a configuration of the nucleons in the single-particle state.
	\author SUN Yazhou
	\date Created: 2020/01/31
	\date Last modified: 2020/03/04 by Sun Yazhou
	\copyright 2020 SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <iomanip>
#include "TAManyBodySD.h"
#include "TAException.h"
#include "TASingleParticleState.h"
#include "TASingleParticleStateManager.h"

using std::cout;
using std::setw;

int TAManyBodySD::kNParticle = -1;

TAManyBodySD::TAManyBodySD(int index, int nParticle, int *SPState)
  : fIndex(index), f2M(0), fEnergy(0.), fBit(){
  kNParticle = nParticle;
  fSPStateArr = new int[kNParticle];
  for(int i = 0; i < kNParticle; i++){
    if(SPState + i) fSPStateArr[i] = SPState[i];
    else TAException::Error("TAManyBodySD",
      "constructor: SPState[%d] is nullptr", i);
  }
  for(int i = 0; i < kNParticle; i++){
    TASingleParticleState *sp = (*this)[fSPStateArr[i]];
    f2M += sp->Get2Mj();
    fEnergy += sp->GetEnergy();
  }
  // conver fSPStateArr to bit into fBit
  fBit.SetBit(fSPStateArr, kNParticle);
} // end of the constructor

TAManyBodySD::~TAManyBodySD(){
  if(fSPStateArr){
    delete [] fSPStateArr;
    fSPStateArr = nullptr;
  }
} // end of the destructor

TASingleParticleState *TAManyBodySD::operator[](int i){
  vector<TASingleParticleState *> spv =
    TASingleParticleStateManager::Instance()->GetSPStateVec();
  if(i >= int(spv.size())){
    TAException::Error("TAManyBodySD",
      "operator[%d]: Subscript of the requested element is out of range.", i);
  }

  TASingleParticleState *sp = spv[i];
  if(!sp){
    TAException::Error("TAManyBodySD",
      "operator[%d]: Required pointer is null.", i);
  }
  return sp;
} // end of member function operator[]

void TAManyBodySD::UpdateSPStateArr(){
  fBit.ConvertToInt(fSPStateArr);
} // end of member function UpdateSPSArr

void TAManyBodySD::UpdateBit(){
  fBit.SetBit(fSPStateArr, kNParticle);
} // end of member function UpdateBit


// self-display
void TAManyBodySD::Print() const{
  cout << std::right;
  cout << "ManyBodySD index: " << setw(4) << fIndex << "   Arr: ";
  for(int i = 0; i < kNParticle; i++) cout << fSPStateArr[i] + 1 << " ";

  cout << "   2M: " << setw(2) << f2M;
  cout << "  energy:" << setw(5) << fEnergy << std::endl;
  cout << std::left;
} // end of member function Print

/// Print the many-body state in bit mode
void TAManyBodySD::PrintInBit() const{
  fBit.PrintInBit();
}

/// \retval fSPStateArr, update with fBit if necessary
int *TAManyBodySD::IntArr() const{
  return fSPStateArr;
}

int TAManyBodySD::GetNParticle(){
  if(kNParticle < 0){
    TAException::Error("TAManyBodySD",
      "GetNParticle: kNParticle: %d is minus, not assigned yet?", kNParticle);
  }
  return kNParticle;
} // end of member function GetNParticle
