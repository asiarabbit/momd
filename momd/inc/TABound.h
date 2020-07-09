/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TABound.h
  \class TABound
  \brief To calculate bound state radial wavefunction of a valence nucleon in
  central potential (Vc+VN+VLS) of a nucleus.    This is a standalone class. It
  is supposed to solicit potential inputs and valence nucleon state quanta from
  users and output the results in arrays to text files.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TABound_h_
#define _TABound_h_

#include <string>

using std::string;

class TABound{
public:
  /// inFile: take user input for potentials and single-particle state
  /// of the valence nucleon
  TABound(const string &inFile = "");
  virtual ~TABound(){}

  /// calculate and output the bound state radial wavefunction in text file
  void ConstructPotential();
  void Bound(double *Rl);

  int Getl() const{ reuturn fl; }

  /// length of the reulting radial wavefunction array
  static const int kNRl = 200;

private:
  /// the Coulomb potential
  double fRC; ///< charge radius
  double fZc, fZv; ///< Z of the core and the valence nucleon
  double fAc, fAv; ///< mass number of the core and the valence nucleon
  /// the nuclear potential: WS form
  double fV0, fR0, fA0; ///< central NN force
  double fVSO, fRSO, fASO; ///< central NN force
  static const int kNV = 200; ///< length of the total potential array
  double *fV; ///< the total potential (VC+VNN+VLS+Vcentrifugal)
  /// single-particle state for the valence nucleon
  int fn, fl, f2j;
};

#endif
