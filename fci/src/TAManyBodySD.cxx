/**
  SUNNY Project, Anyang Normal University, IMP-CAS
	\file TAManyBodySD.cxx
	\class TAManyBodySD
	\brief Slater determinant (SD) class for many-body problems. Each SD represents
	a configuration of the nucleons in the single-particle state.
	\author SUN Yazhou
	\date Created: 2020/01/31
	\date Last modified: 2021/02/10 by Sun Yazhou
	\copyright 2020 SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include "TAManyBodySD.h"
#include "TAException.h"
#include "TASPState.h"
#include "TASPStateManager.h"

using std::cout;
using std::setw;

TAManyBodySD::TAManyBodySD(unsigned long bit) : fBit(bit){
} // end of the constructor

TAManyBodySD::~TAManyBodySD(){} // end of the destructor


void TAManyBodySD::GetSPStateArr(int *p) const{
  const int n = TASPStateManager::Instance()->GetSPStateVec().size();
  int j = 0;
  for(int i = 0; i < n; i++) if(fBit.test(i)) p[j++] = i;
}

short TAManyBodySD::Get2M() const{
  short total2M = 0;
  int n = fBit.count();
  for(int i = n; i--;) total2M += (*this)[i]->Get2Mj();
  return total2M;
}
double TAManyBodySD::GetEnergy() const{
  double e = 0.;
  int n = fBit.count();
  for(int i = n; i--;) e += (*this)[i]->GetEnergy();
  return e;
}

/// note that this method works only if piared states are next to each other,
/// and SPStates are ordered in ManyBodySD
bool TAManyBodySD::IsPaired() const{
  // obtain the sp state array
  const int n = fBit.count();
  int p[n]{};
  GetSPStateArr(p);

  // identify broken pairs //
  int i = 0, nSingle = 0;
  while(i < n){
    if(i < n-1 && !(p[i]%2) && p[i+1]-p[i] == 1) i += 2; // p[i] even and p[i+1] odd
    else{
      nSingle++;
      i++;
    }
  } // end while
  if(nSingle <= 1) return true;
  return false;
} // end of member function IsPaired

/// \retval sp state of the i-th particle
TASPState *TAManyBodySD::operator[](int i) const{
  static vector<TASPState *> &spv =
    TASPStateManager::Instance()->GetSPStateVec();
  int np = -1, j, n = spv.size();
  if(i >= n) TAException::Error("TAManyBodySD",
    "operator[%d]: Subscript of the requested element is out of range.", i);
  // obtain the sp state array /
  for(j = 0; j < n; j++){
    if(fBit.test(j)) np++;
    if(np == i) break;
  }
  TASPState *sp = spv[j];
  if(!sp) TAException::Error("TAManyBodySD", "operator[%d]: Required pointer is null.", i);
  return sp;
} // end of member function operator[]

// self-display
void TAManyBodySD::Print() const{
  const int n = fBit.count();
  int p[n]{};
  GetSPStateArr(p);

  cout << std::right;
  cout << "ManyBodySD Arr: ";
  for(int i = 0; i < n; i++) cout << p[i] + 1 << " ";

  cout << "   2M: " << setw(2) << Get2M();
  cout << "  energy:" << setw(5) << GetEnergy() << std::endl;
  cout << std::left;
} // end of member function Print

/// Print the many-body state in bit mode
void TAManyBodySD::PrintInBit() const{
  fBit.PrintInBit();
}
