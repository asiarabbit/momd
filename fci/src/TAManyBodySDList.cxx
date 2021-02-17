/**
	\file TAManyBodySDList.C
	\class TAManyBodySDList
	\brief A list to store many-body Slater determinants, classified by M, which is
	the 3rd component jz of the total angular momentum. So M is the same for each
	member of this list.
	\author SUN Yazhou
	\date Created: 2020/01/31
	\date Last modified: 2021/02/06 by SUN Yazhou
	\copyright 2020 SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <algorithm>
#include "TAManyBodySDList.h"
#include "TAManyBodySD.h"
#include "TAException.h"
#include "TABit.h"
#include "TASPStateManager.h"

using std::cout;
using std::endl;

TAManyBodySDList::TAManyBodySDList(short twoM) : f2M(twoM){
  fManyBodySDTree = new TAMBSDTree("mbsdList", "Many-Body Slater Determinant Tree");
}

TAManyBodySDList::~TAManyBodySDList(){
  if(fManyBodySDTree){
    delete fManyBodySDTree; fManyBodySDTree = nullptr;
  }
} // end of the destructor

bool TAManyBodySDList::MBSDFilter(unsigned long long bit) const{
  TAManyBodySD m(bit);
  return
    // m.IsPaired() &&
    m.Get2M() == f2M;
}

void TAManyBodySDList::Add(unsigned long long bit){
  if(MBSDFilter(bit)) fManyBodySDTree->Fill(bit);
}

void TAManyBodySDList::Add(int *spsArr, int nParticle){
  TABit bit; bit.SetBit(spsArr, nParticle);
  Add(bit.to_ulong());
}

unsigned long long TAManyBodySDList::Bit(unsigned long long mbsdIndex){
  return fManyBodySDTree->GetMBSDInBit(mbsdIndex);
}

void TAManyBodySDList::Print(){
  cout << "Print many-body basis set for M-scheme where 2M = ";
  cout << f2M << endl;
  int n = GetNBasis();
  if(n > 30) n = 30;
  for(int i = 0; i < n; i++) TAManyBodySD((*this)[i]).Print();
  cout << "Totally there're " << n;
  cout << " many-body Slater determinants in the list with f2M = ";
  cout << f2M << "." << endl;
}

/// Print all the mbsd-s in bit mode
void TAManyBodySDList::PrintInBit(){
  cout << "Print many-body basis set for M-scheme where 2M = ";
  cout << f2M << " in bit mode" << endl;
  int n = GetNBasis();
  if(n > 30) n = 30;
  for(int i = 0; i < n; i++) TAManyBodySD((*this)[i]).PrintInBit();
  cout << "Totally there're " << n;
  cout << " many-body Slater determinants in the list with f2M = ";
  cout << f2M << "." << endl;
}

unsigned long long TAManyBodySDList::GetNBasis() const{
  return fManyBodySDTree->GetEntries();
}

///< note that index is the mbsd index
TASPState *TAManyBodySDList::GetSPState(unsigned long long mbsdIndex, int i){
  static vector<TASPState *> &spv =
    TASPStateManager::Instance()->GetSPStateVec();
  if(i >= int(spv.size()))
    TAException::Error("TAManyBodySDList",
      "GetSPState: Subscript of the requested sp orbit is out of range.", i);

  TABit b(fManyBodySDTree->GetMBSDInBit(mbsdIndex));
  int np = -1, j, n = b.size();
  for(j = 0; j < n; j++){
    if(b.test(j)) np++;
    if(np == i) break;
  }
  TASPState *sp = spv[j];
  if(!sp) TAException::Error("TAManyBodySDList", "GetSPState: Required sp is null.");
  return sp;
} // end of member function GetSPState

/// \retval <rr|a+_p * a_q|cc>
int TAManyBodySDList::Integral(int rr, int p, int q, int cc){
  TABit rBit((*this)[rr]);
  rBit.Annhilate(p); if(!rBit.GetPhase()) return 0.;

  TABit cBit((*this)[cc]);
  cBit.Annhilate(q); if(!cBit.GetPhase()) return 0.;

  return rBit*cBit;
} // end of member function Integral(rr,p,q,cc);

/// \retval <rr|a+_p*a+_q * a_s*a_r|cc>
int TAManyBodySDList::Integral(int rr, int p, int q, int s, int r, int cc){
  TABit rBit((*this)[rr]);
  rBit.Annhilate(p); if(!rBit.GetPhase()) return 0.;
  rBit.Annhilate(q); if(!rBit.GetPhase()) return 0.;

  TABit cBit((*this)[cc]);
  cBit.Annhilate(r); if(!cBit.GetPhase()) return 0.;
  cBit.Annhilate(s); if(!cBit.GetPhase()) return 0.;

  return rBit*cBit;
} // end of member function Integral(rr,p,q,r,s,cc);
/// \retval <rr|a+_p*a+_q*a+_r * a_u*a_t*a_s|cc>
int TAManyBodySDList::Integral(int rr, int p, int q, int r, int u, int t, int s, int cc){
  TABit rBit((*this)[rr]);
  rBit.Annhilate(p); if(!rBit.GetPhase()) return 0.;
  rBit.Annhilate(q); if(!rBit.GetPhase()) return 0.;
  rBit.Annhilate(r); if(!rBit.GetPhase()) return 0.;

  TABit cBit((*this)[cc]);
  cBit.Annhilate(s); if(!cBit.GetPhase()) return 0.;
  cBit.Annhilate(t); if(!cBit.GetPhase()) return 0.;
  cBit.Annhilate(u); if(!cBit.GetPhase()) return 0.;

  return rBit*cBit;
} // end of member function Integral(rr,p,q,r,s,t,u,cc);
