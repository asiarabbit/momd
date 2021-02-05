/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TABit.h
  \class TABit
  \brief many-body Slater determinant in bit representation. This class is coined
  to be a member of a many-body state. Also incorporated are create and annhilate
  operations. This is also partly used as a light-weight MBSD. Inherited from bitset.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/11
  \date Last modified: 2021/02/04 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TABit_h_
#define _TABit_h_

#include <bitset>

using std::bitset;
static const int NBIT = 32;

class TABit : public bitset<NBIT>{
public:
  TABit(ulong n = 0); ///< n: bit arragement in decimal format
  TABit(const TABit &bit); ///< copy constructor
  TABit &operator=(const TABit &bit); ///< assignment constructor
  /// \retval: <*this|bit> with phase updated
  int operator*(const TABit &bit) const;
  virtual ~TABit();


  /// \brief conversion from int array to bit, e.g. {0,3,7} -> 10010001000...
  /// \param np: number of particles, length of arr
  void SetBit(const int *arr, int np);
  void Reset(); ///< \brief set the bit to zero
  /// \brief conversion from bit to int array, the inverse of SetBit method
  void ConvertToInt(int *arr);
  /// apply a+_p operator and alter fBit and fPhase accordingly \retval *this
  TABit &Create(int p);
  /// apply a_p operator and alter fBit and fPhase accordingly \retval *this
  TABit &Annhilate(int p);

  short GetPhase() const{ return fPhase; }
  void PrintInBit() const; ///< print in bit

private:
  short fPhase; ///< 0, +-1
};

#endif
