/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFNN.h
  \class TAFNN
  \brief This class represents nucleon-nucleon scattering amplitude, extracted
  from fitting of pp, pn and nn scattering data over a wide range of energies.
  fNN(q)=kNN/4pi*sigmaNN(i+alphaNN)*exp(-betaNN*q^2).
  This is a static class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/25 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAFNN_h_
#define _TAFNN_h_

class TAFNN{
public:
  enum ParOption{kHoriuchi = 0, kLenziRay};
  TAFNN();
  virtual ~TAFNN();

  /// isospin-averaged nucleon-nucleon scattering cross section
  /// for lab energy Ek.
  /// \param isPauli 0: free nn scattering; 1: Pauli-blocking corrected
  /// P: projectile, T: target. a, z: A(atomic number) and Z(nuclear charge)
  /// \retval in fm^2
  static double GetSigmaNN(double pEk, int zP, int nP, int zT, int nT,
    bool isPauli = false);
  /// fNN(q)=kNN/4pi*sigmaNN(i+alphaNN)*exp(-betaNN*q^2)
  /// alphaNN at arbitrary Ek is interpolated from array alphaNN
  static double GetAlphaNN(double pEk, int zP, int nP, int zT, int nT);
  /// betaNN at arbitrary Ek is interpolated from array alphaNN
  static double GetBetaNN(double pEk, int zP, int nP, int zT, int nT);
  static void SetParOption(ParOption opt){ kParOpt = opt; }

  static ParOption kParOpt; ///< parmeterization option for alphaNN and betaNN

protected:
  /// P: projectile, T: target
  /// \retval (vpp*(nP*zP+nT*zT) + vnp*(nP*zP+nT*zT)) / (aP*aT)
  static double IsospinAverage(double vpp, double vpn, int zP, int nP,
      int zT, int nT);
  /// Horiuchi, et al., PRC.75.044607. E: [30, 1000], N=16
  static double AlphaBetaHoriuchi(double pEk, bool isAlpha = true);
  /// Ray, et al., PRC.20.1857. E: [100, 2200], N=10
  /// Lenzi, et al., PRC.40.2114. E: [10, 94], N=7
  static double AlphaBetaLenziRay(double pEk, int zP, int nP, int zT, int nT,
      bool isAlpha = true);
};

#endif
