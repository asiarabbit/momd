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

using std::cout;
using std::endl;
using std::remove_if;

TAManyBodySDList::TAManyBodySDList(short twoM) : f2M(twoM){
}

TAManyBodySDList::~TAManyBodySDList(){}

void TAManyBodySDList::Add(TAManyBodySD *mbsd){
  if(!mbsd) TAException::Error("TAManyBodySD", "Add: Input is a nullptr");
  fManyBodySDVec.push_back(mbsd);
  //mbsd->SetIndex(fManyBodySDVec.size() - 1);
}

/// remove those with broken pairs, i.e. >2 single particles
void TAManyBodySDList::Pairing(){
  auto &p = fManyBodySDVec;
  p.erase(remove_if(p.begin(), p.end(), [](TAManyBodySD *l){ return !(l->IsPaired()); }), p.end());
} // end member function Pairing

void TAManyBodySDList::Print() const{
  cout << "Print many-body basis set for M-scheme where 2M = ";
  cout << f2M << endl;
  for(TAManyBodySD *mp : fManyBodySDVec) mp->Print();
  cout << "Totally there're " << fManyBodySDVec.size();
  cout << " many-body Slater determinants in the list with f2M = ";
  cout << f2M << "." << endl;
}

/// Print all the mbsd-s in bit mode
void TAManyBodySDList::PrintInBit() const{
  cout << "Print many-body basis set for M-scheme where 2M = ";
  cout << f2M << " in bit mode" << endl;
  for(TAManyBodySD *mp : fManyBodySDVec) mp->PrintInBit();
  cout << "Totally there're " << fManyBodySDVec.size();
  cout << " many-body Slater determinants in the list with f2M = ";
  cout << f2M << "." << endl;
}

TAManyBodySD *TAManyBodySDList::operator[](int i) const{
  if(i > int(fManyBodySDVec.size()) || i < 0){
    TAException::Error("TAManyBodySD",
      "operator[]: index i: %d out of range.", i);
  }
  return fManyBodySDVec[i];
} // end of member function operator[]

/// \retval <rr|a+_p * a_q|cc>
int TAManyBodySDList::Integral(int rr, int p, int q, int cc) const{
  TABit rBit = (*this)[rr]->Bit();
  rBit.Annhilate(p); if(!rBit.GetPhase()) return 0.;

  TABit cBit = (*this)[cc]->Bit();
  cBit.Annhilate(q); if(!cBit.GetPhase()) return 0.;

  return rBit*cBit;
} // end of member function Integral(rr,p,q,cc);

/// \retval <rr|a+_p*a+_q * a_s*a_r|cc>
int TAManyBodySDList::Integral(int rr, int p, int q, int s, int r, int cc) const{
  TABit rBit = (*this)[rr]->Bit();
  rBit.Annhilate(p); if(!rBit.GetPhase()) return 0.;
  rBit.Annhilate(q); if(!rBit.GetPhase()) return 0.;

  TABit cBit = (*this)[cc]->Bit();
  cBit.Annhilate(r); if(!cBit.GetPhase()) return 0.;
  cBit.Annhilate(s); if(!cBit.GetPhase()) return 0.;

  return rBit*cBit;
} // end of member function Integral(rr,p,q,r,s,cc);
/// \retval <rr|a+_p*a+_q*a+_r * a_u*a_t*a_s|cc>
int TAManyBodySDList::Integral(int rr, int p, int q, int r, int u, int t,
    int s, int cc) const{
  TABit rBit = (*this)[rr]->Bit();
  rBit.Annhilate(p); if(!rBit.GetPhase()) return 0.;
  rBit.Annhilate(q); if(!rBit.GetPhase()) return 0.;
  rBit.Annhilate(r); if(!rBit.GetPhase()) return 0.;

  TABit cBit = (*this)[cc]->Bit();
  cBit.Annhilate(s); if(!cBit.GetPhase()) return 0.;
  cBit.Annhilate(t); if(!cBit.GetPhase()) return 0.;
  cBit.Annhilate(u); if(!cBit.GetPhase()) return 0.;

  return rBit*cBit;
} // end of member function Integral(rr,p,q,r,s,t,u,cc);
