/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TABit.cxx
  \class TABit
  \brief 128-bit for many-body basis in bit representation. This class is coined
  to be a member of a many-body state. Also incorporated are create and
  annhilate operations. This is also partly used as a light-weight MBSD.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/11
  \date Last modified: 2020/02/11 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <cstring>
#include <iostream>
#include <bitset> // for printing in bit
#include "TABit.h"
#include "TAException.h"

TABit::TABit() : fBit{0}, fPhase(1){}

/// copy constructor
TABit::TABit(const TABit &bit){
  *this = bit;
} // end of the copy constructor

/// assignment constructor
TABit &TABit::operator=(const TABit &bit){
  if(&bit == this) return *this;
  memcpy(fBit, bit.fBit, sizeof(fBit));
  fPhase = bit.fPhase;
  return *this;
} // end of the assignment constructor

TABit::~TABit(){}

/// \brief conversion from int array to bit, e.g. {0,3,7}: 10010001000...
/// \param np: length of arr; bit[8]: 0-8: ascending order
void TABit::SetBit(const int *arr, int np){
  Reset(); // set fBit array to zero
  for(int i = 0; i < np; i++){
    if(arr[i] < 0 || arr[i] >= int(sizeof(fBit)*8)){
      TAException::Error("TABit",
        "SetBit: The input SP state ouf of range, arr[%d]: %d", i, arr[i]);
    }
    fBit[arr[i]/32] += 1 << (arr[i] % 32);
  } // end for over i
} // end member function Bit

/// \brief set the bit to zero
void TABit::Reset(){
  memset(fBit, 0, sizeof(fBit));
  fPhase = 1;
}

/// \brief conversion from bit to int array
void TABit::ConvertToInt(int *arr){
  static const int nbit = sizeof(fBit) * 8;
  int np = 0; // number of particles
  for(int i = 0; i < nbit; i++){
    if(fBit[i/32] & (1 << i%32)){ // single-particle state i occupied
      try{
        // \NOTE that mbsd-s are automatically put in ascending order in arr
        arr[np++] = i;
      }
      catch(...){
        TAException::Error("TABit",
          "ConvertToInt: Illegal memory access, np: %d", np);
      }
    } // end if
  } // end for over i
} // end member function ConvertToInt

/// apply a+_p operator and alter fBit and fPhase accordingly \retval *this
/// if the result of operation is zero, zero will be stored in phase
/// instead of the bit array
TABit &TABit::Create(int p){
  if(!fPhase) return *this; // zero state, waste of time

  static const int nbit = sizeof(fBit) * 8;
  if(p >= nbit || p < 0){
    TAException::Error("TABit",
      "Create: SP state %d to create a particle on \
is out of range. nbit: %d", p, nbit);
  }
  int np = 0; // number of particles
  for(int i = 0; i < p; i++){
    if(fBit[i/32] & (1 << i%32)) np++; // single-particle state i occupied
  } // end for over i
  if(fBit[p/32] & (1 << p%32)){ // single-particle state p is occupied
    fPhase = 0; // occupied in SPS p, cannot create on an occupied SPS
    return *this; // Pauli's exclusion principle
  }
  else fBit[p/32] += (1 << p%32); // a particle on SPS p created
  // assign the phase: phase*(-)^np
  if(np%2) fPhase *= -1;
  return *this;
} // end of member function Create

/// apply a_p operator and alter fBit and fPhase accordingly \retval *this
TABit &TABit::Annhilate(int p){
  if(!fPhase) return *this; // zero state, waste of time

  static const int nbit = sizeof(fBit) * 8;
  if(p >= nbit || p < 0){
    TAException::Error("TABit",
      "Annhilate: SP state %d to annhilate a particle on \
is out of range. nbit: %d", p, nbit);
  }
  int np = 0; // number of particles
  for(int i = 0; i < p; i++){
    if(fBit[i/32] & (1 << i%32)) np++; // single-particle state i occupied
  } // end for over i
  if(fBit[p/32] & (1 << p%32)){ // single-particle state p is occupied
    fBit[p/32] -= (1 << p%32); // a particle on SPS p annhilated
  }
  else{
    fPhase = 0; // empty in SPS p, cannot annhilate on an empty SPS
    return *this; // Pauli's exclusion principle
  }
  // assign the phase: phase*(-)^np
  if(np%2) fPhase *= -1;
  return *this;
} // end of member function Annhilate

/// \retval: <*this|bit> phase included
int TABit::operator*(const TABit &bit) const{
  // to tell if either is zero phase
  int phase = fPhase * bit.fPhase;
  if(!phase) return 0;

  // to tell if the two states are orthogonal or parallel
  static const int nword = sizeof(fBit) / sizeof(unsigned);
  for(int i = 0; i < nword; i++){
    try{
      if(fBit[i] != bit.fBit[i]) return 0; // orthogonal
    }
    catch(...){
      TAException::Error("TABit",
        "operator*: Illegal memory access, i: %d", i);
    }
  } // end for over words of fBit
  return phase; // the two states are parallel
} // end of member function operator*

/// Print in bit
void TABit::PrintInBit() const{
  static const int nbit = sizeof(fBit) * 8; // 1 byte = 8 bits
  std::bitset<nbit> b(0);
  for(int i = 0; i < nbit; i++){
    if(fBit[i/32] & (1 << i%32)) b.set(i);
  }
  std::cout << b << std::endl;
} // end of function Print
