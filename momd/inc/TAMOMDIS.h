/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMOMDIS.h
  \class TAMOMDIS
  \brief This is a global class, to take control of the whole flow of the
  program. It is also responsible for generating various momentum distributions.
  This is supposed to be a singleton class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <string>

using std::string;

#ifndef _TAMOMDIS_h_
#define _TAMOMDIS_h_

class TAMOMIS{
public:
  static TAMOMDIS *Instance();
  virtual ~TAMOMDIS();

  void Go();
  void Prepare(); ///< solve Rl, generate and fit S-matrix
  void Parallel(); ///< calculate dsigma/dkz

protected:
  TAMOMDIS(); ///< the constructor

  static TAMOMDIS *kInstance;
  /// specify the config file folder
  string fConfigDir;
  double *fMOMPara; ///< parallel momentum (stripping) distribution
  double fSigmaStr; ///< stripping c.s.
  double fSigmaDiff; ///< diffraction c.s.
  int fl; ///< angular momentum of valence nucleon sp state

  /// the component to calculate momtum distribution of the core
  TABound *fBound;
  double *fRl; ///< the radial wavefunction
  TASMatrix *fSMatrix;
  /// the fitting result of the S-matrix Gaussian expansion
  TAComplex *fAlphaj; ///< \sum_j(  alphaj*exp(-b^2/betaj^2)  )
  double *fBetaj;
};

#endif
