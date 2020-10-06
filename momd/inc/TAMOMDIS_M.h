/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMOMDIS_M.h
  \class TAMOMDIS_M
  \brief Calculate residue momentum distribution with angular momentum component
  m specified. This is a class to assist class TAMOMDIS.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/10/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMOMDIS_M_h_
#define _TAMOMDIS_M_h_

#include <complex>
#include "TAMOMDIS.h"

typedef std::complex<double> cdouble;

class TAMOMDIS;

class TAMOMDIS_M{
public:
  TAMOMDIS_M(TAMOMDIS *mom);
  virtual ~TAMOMDIS_M();

  /// returning c.s. for m, mom is assigned to momStr
  double ParallelStr(int l, int m, double *momStr);
  void ParallelStr(int l, int m, double kz, double &momStr); ///<\retval dsig_str/dkz
  // calculate the same result as ParallelStr, but with a somewhat more direct method //
  // i.e. sc is integrated over phi without using sum of gaussians to approximate sc
  void ParallelStr1(int l, int m, double kz, double &momStr);
  /// the same as ParallelStr, but for diffraction dissociation
  double ParallelDiff(int l, int m, double *momDiff);
  void ParallelDiff(int l, int m, double kz, double &momDiff); ///<\retval dsig_diff/dkz

  double *GetMomArr() const{ return fMomArr; }

  static const int kNmom = 200; ///< number of points in momentum distribution

protected:
  double GetRl(double r); ///< the radial wavefunction
  cdouble GetSc(double b){ return fMOM->GetSc(b); }
  cdouble GetSn(double b){ return fMOM->GetSn(b); }
  // integrate sc over phi, the azimuthal angle of r w.r.t. b
  double IntToRhoBStr(double rho, double b, double kz, int l, int m);

  TAMOMDIS *fMOM;
  double *fMomArr; ///< momentum array for the momentum distribution
  // parameters all from fMOM->fSc
  int fNG; // number of Guassians in Gaussian fit of the SMatrix
  double *fAlphajR, *fAlphajI, *fBetaM2, fRL; // fBetaM2 = 1/betaj^2
};

#endif
