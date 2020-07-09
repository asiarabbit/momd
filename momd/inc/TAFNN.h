/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFNN.h
  \class TAFNN
  \brief This class represents nucleon-nucleon scattering amplitude, extracted
  from fitting of pp, pn and nn scattering data over a wide range of energies.
  fNN(q)=kNN/4pi*sigmaNN(i+alphaNN)*exp(-betaNN*q^2)
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAFNN_h_
#define _TAFNN_h_

class TAComplex;

class TAFNN{
public:
  TAFNN();
  virtual ~TAFNN();

  /// \param kNN: wave number of n-n scattering
  TAComplex GetFNN(double kNN);

private:
  static const int kNE = 10; ///< number of c.s. samplings w.r.t. energy
  double fEnergy[kNE];
  /// n-n scattering cross section
  double fSigmapp[kNE];
  double fSigmapn[kNE];
  // alphaNN and betaNN
  double fAlphapp[kNE];
  double fAlphapn[kNE];
  double fBetapp[kNE];
  double fBetapn[kNE];
};
