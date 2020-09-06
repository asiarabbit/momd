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
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TABound_h_
#define _TABound_h_

#include <string>
#include <vector>

using std::string;
using std::vector;

class TABound{
public:
  /// inFile: take user input for potentials and single-particle state
  /// of the valence nucleon
  TABound(const vector<double> &vv);
  virtual ~TABound();

  /// calculate the bound state radial wavefunction in fV
  void ConstructPotential();
  void Bound(); ///< calculate and output the bound state radial wavefunction in fRl
  const double *GetRl(); ///< output the solved wavefunction

  int Getl() const{ return fl; }

  static const int kNRl = 200; /// length of the reulting radial wavefunction array
  static const int kNV = 200; ///< length of the total potential array

private:
  double fZc, fAc; ///< identity of the core
  double fZv, fAv; ///< identity of the valence nucleon
  /// single-particle state for the valence nucleon
  int fn, fl, f2j;
  /// the nuclear potential: WS form
  double fV0, fR0, fA0; ///< central NN force
  double fVS, fRS, fAS; ///< central NN force
  double fRC; ///< charge radius for Coulomb potential

  double *fV; ///< the total potential (VC+VNN+VLS+Vcentrifugal)
  double *fRl; ///< the resulting radial wavefunction
};

#endif
