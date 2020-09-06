/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASMatrix.cxx
  \class TASMatrix
  \brief A collection of methods to fit data to certain models.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/27
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cmath>
#include "TAFit.h"
#include "TAMath.h"
#include "TAException.h"

/// Given a set of data x[0...ndata-1] and y[0...ndata-1] with standard deviations
/// sigy[0...ndata-1]. Fit them to a straight line by minimizing chi2. Returned
/// are a and b, and their respective probable uncertainties.siga and sigb, the
/// chi-square chi2 and the goodness-of-fit probability q (that the fit would
/// have chi2 this large or larger under the model given by the fitting). If
/// dy[0] <= 0., then the standard deviations are assumed to be unavailable: q
/// is returned as 1.0 and the normalization of chi2 is to unit standard deviation
/// on all points.
/// Ref. Numerical Recipes in C, p689.
void TAFit::LinearFit(const double *x, const double *y, const double *dy, int ndata,
    double &a, double &b, double &siga, double &sigb, double &chi2, double &q){
  // stt=\sum_i{1/sigi*(xi-sx/s), sx=\sum_i{xi/sigi^2}, sy=\sum_i{yi/sigi^2}
  double s = 0., sx = 0., sy = 0., stt = 0.;
  bool isSigYKnown = dy[0] > 0.;
  // accumulate sums ...
  if(isSigYKnown){
    s = 0.;
    for(int i = 0; i < ndata; i++){ // ...with weights
      double wt = 1. / (dy[i]*dy[i]); // weights
      s  += wt;
      sx += wt*x[i];
      sy += wt*y[i];
    } // end for
  } // end if
  else{ // ...or without weights
    for(int i = 0; i < ndata; i++){
      sx += x[i];
      sy += y[i];
    } // end for over i
    s = ndata;
  } // end else
  double sxos = sx / s; // sx over s
  // calculate the slope b
  b = 0.;
  if(isSigYKnown){
    for(int i = 0; i < ndata; i++){
      double t = (x[i] - sxos)/dy[i];
      stt += t*t;
      b += t*y[i]/dy[i];
    } // end for
  } // end if
  else{
    for(int i = 0; i < ndata; i++){
      double t = x[i] - sxos;
      stt += t*t;
      b += t*y[i];
    } // end for over i
  } // end else

  // solve for a, b, siga and sigb
  b /= stt;
  a = (sy-sx*b)/s;
  siga = sqrt((1.+sx*sx/(s*stt))/s);
  sigb = sqrt(1./stt);

  // calculte chi2
  chi2 = 0.; q = 0.;
  double nu = ndata - 2; // degree of freedom for chi2 distribution
  if(isSigYKnown){
    for(int i = 0; i < ndata; i++) chi2 += TAMath::sqr((y[i]-a-b*x[i])/dy[i]);
    q = TAMath::gammaq(0.5*nu, 0.5*chi2);
  } // end if
  // for unweighted data, evaluate typical sig (which is assumed unit before)
  // using chi2, then update siga and sigb
  else{
    for(int i = 0; i < ndata; i++) chi2 += TAMath::sqr(y[i]-a-b*x[i]);
    // chi2~nu, chi2*sig^2=\sum{{yi-y(xi)}^2}~nu*sig^2
    double sigData = sqrt(chi2 / nu);
    siga *= sigData;
    sigb *= sigData;
  } // end else
} // end of member function LinearFit

/// General Least Square Method fitting - Given a set of data x[0...ndata-1]
/// and y[0...ndata-1] with standard deviations dy[0...ndata-1]. Use chi2
/// minimization to fit for some or all of the coefficients a[0...ma-1] of a
/// function taht depends linearly on a: y =\sum_i{ai*func_i(x)}. The input
/// boolean array isFit[0...ma-1] indicates by false those components of a that
/// should be held fixed at their input values. The program returns values for
/// a[0...ma-1], chi2, and the covariance matrix covar[0...ma-1][0...ma-1].
/// Parameters held fixed will return zero covariances. Users should supply a
/// routine fun(x, funci, ma) that returns the ma basis functions evaluated at
/// x in the array funci[0...ma-1]
/// Ref. Numerical Recipes in C: p675
void TAFit::LSMFit(const double *x, const double *y, const double *dy, int ndata,
    int ma, double *a, bool *isFit, matrix &covar, double &chi2,
    void (*fun)(double x, double *funci, int ma, const double *p), const double *p){
  int nFit = 0; // number of parameters to be fitted
  for(int i = 0; i < ma; i++) if(isFit[i]) nFit++;
  if(!nFit) TAException::Error("TAFit", "LSMFit: No parameters to be fitted.");
  double beta[nFit]{}, funci[ma]{}; // beta=A^T.b, funci=func_i(x)
  for(int i = 0 ; i < nFit; i++) for(int j = 0 ; j < nFit; j++) covar[i][j] = 0.;

  // loop over data to accumulate coefficients of the normal equations
  int jj = 0, kk = 0; // jj, kk: subscripts for beta and covar
  for(int i = 0; i < ndata; i++){
    fun(x[i], funci, ma, p); // assign array funci
    double yi = y[i];
    // downgrade the problem to nFit dimensions
    if(nFit < ma) for(int j = 0; j < ma; j++) if(!isFit[j]) yi -= a[j]*funci[j];
    const double wt0 = TAMath::sqr(dy[i]);
    for(int j = jj = 0; j < ma; j++){
      if(isFit[j]){
        const double wt = funci[j]/wt0;
        for(int k = kk = 0; k < j; k++) if(isFit[k]) covar[jj][kk++] += funci[k]*wt;
        beta[jj++] += yi*wt;
      } // end outer if
    } // end for over j
  } // end loop over data points
  // fill in covar above the diagonal from symmetry
  for(int i = 1; i < nFit; i++) for(int j = i + 1; j < nFit; j++) covar[i][j] = covar[j][i];
  // solve matrix equation covar*a=b using Gauss-Jordan elimination
  // covar would be inversed inplace, and the solution would be stored in beta
  TAMath::GaussJordan(covar, nFit, beta);
  // assign the solution to array a
  for(int i = jj = 0; i < ma; i++) if(isFit[i]) a[i] = beta[jj++];
  // calculate chi2
  for(int i = 0, chi2 = 0.; i < ndata; i++){
    fun(x[i], funci, ma, p);
    double sum = 0.;
    for(int j = 0; j < ma; j++) sum += a[j]*funci[j];
    chi2 += TAMath::sqr((y[i]-sum)/dy[i]);
  } // end for over i
  // expand covar to include the fixed parameters and in a[0...ma-1] order
  TAMath::FillCovar(covar, ma, isFit, nFit);
} // end of member function LSMFit
