/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASMatrix.h
  \class TASMatrix
  \brief A collection of methods to fit data to certain models.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/27
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAFit_h_
#define _TAFit_h_

#include "TAMatrix.h"

class TAFit{
public:
  TAFit(){}
  virtual ~TAFit(){}

  /// Given a set of data x[0...ndata-1] and y[0...ndata-1] with standard deviations
  /// sigy[0...ndata-1]. Fit them to a straight line by minimizing chi2. Returned
  /// are a and b, and their respective probable uncertainties.siga and sigb, the
  /// chi-square chi2 and the goodness-of-fit probability q (that the fit would
  /// have chi2 this large or larger under the model given by the fitting). If
  /// dy[0] < 0., then the standard deviations are assumed to be unavailable: q
  /// is returned as 1.0 and the normalization of chi2 is to unit standard deviation
  /// on all points.
  static void LinearFit(const double *x, const double *y, const double *dy, int ndata,
      double &a, double &b, double &siga, double &sigb, double &chi2, double &q);

  /// General Least Square Method fitting - Given a set of data x[0...ndata-1]
  /// and y[0...ndata-1] with standard deviations dy[0...ndata-1]. Use chi2
  /// minimization to fit for some or all of the coefficients a[0...ma-1] of a
  /// function taht depends linearly on a: y =\sum_i{ai*func_i(x)}. The input
  /// boolean array isFit[0...ma-1] indicates by false those components of a that
  /// should be held fixed at their input values. The program returns values for
  /// a[0...ma-1], chi2, and the covariance matrix covar[0...ma-1][0...ma-1].
  /// Parameters held fixed will return zero covariances. Users should supply a
  /// routine fun(x, funci, ma, p) that returns the ma basis functions evaluated at
  /// x in the array funci[0...ma-1]. p here is for possible parameters in the Gaussians
  /// Ref. Numerical Recipes in C: p675
  static void LSMFit(const double *x, const double *y, const double *dy, int ndata,
      int ma, double *a, bool *isFit, matrix &covar, double &chi2,
      void (*fun)(double x, double *funci, int ma, const double *p), const double *p = nullptr);
};

#endif
