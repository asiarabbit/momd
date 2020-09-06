/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASMatrix.h
  \class TASMatrix
  \brief This class calculates S-matrix exp(i*chi), from given eikonal phase chi.
  Mostly, this class fits the S-matrix with sum of Gaussians, and output the
  fitting result alphaj, i.e., exp(i*chi)=\sum_j{(alphajR+i*alphajI)*exp(-b^2/betaj)},
  where betaj=RL/j.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASMatrix_h_
#define _TASMatrix_h_

#include <complex>
#include <string>
#include "TAException.h"

using cdouble = std::complex<double>;
using std::string;

class TAEikonalPhase;

class TASMatrix{
public:
  // ep, ng, nb, rl are initializations for
  // fEikonalPhase, fNG fNB and fRL respectively
  TASMatrix(int zP, int aP, int zT, int aT, double ek);
  virtual ~TASMatrix();

  /// calculate and output the S-matrix
  cdouble SMatrix(double b) const;
  ///< real and imag part of S-matrix, in an array (0, kBMax) with step h=kBMax/fNB
  void SMatrix(double *b, double *smr, double *smi) const;
  /// fit the S-matrix to an expansion of Gaussians
  /// the results are stored in input arrays alphaj and betaj and members fAlphajR,I
  void GaussianFit();

  double *GetAlphajR();
  double *GetAlphajI();
  void SetRL(double rl){ fRL = rl; }
  /// gaussians used for S-Matrix expansion, \retval array funci[ma] evaluated at x
  static void gaus(double x, double *funci, int ma, const double *p);

protected:

  TAEikonalPhase *fEikonalPhase; ///< contains everything needed to calculate eikonal phase
  int fNB; ///< number of b sampling for eikonal phase
  int fNG; ///< number of Gaussians used in the expansion
  ///< beta[j] = fRL/j; fRL is typical of nuclear size, e.g. RL~50fm for 15C+B9
  double fRL;
  /// the fitting results
  double *fAlphajR, *fAlphajI;
  /// as the name indicates, unit: fm
  double fBMax;
};

#endif
