/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASMatrix.h
  \class TASMatrix
  \brief This class calculate the eikonal S-matrix. It takes nucleus densities
  and nucleon-nucleon scattering amplitude (fNN) as input. The output is stored
  in arrays.   Since the S-matrix is fitted with sum of Gaussians, this class
  takes care of the fitting as well.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASMatrix_h_
#define _TASMatrix_h_

class TAComplex;

class TASMatrix{
public:
  TASMatrix();
  virtual ~TASMatrix();

  /// calculate and output the S-matrix
  void SMatrix(TAComplex *smat = nullptr);
  /// fit the S-matrix to an expansion of Gaussians
  /// the results are stored in member arrays alphaj and betaj
  void GaussianFit(double *alphaj, double *betaj);

  void GetAlphaj(double *alphaj);
  void GetBetaj(double *betaj);

  static const int kNG = 21; ///< number of Gaussians used in the expansion

private:
  static const int kNR = 200; ///< number of r sampling in nuclear density
  /// input nulcear densities and n-n scattering amplitude
  TAComplex *fRhoP; ///< the 2-D Fourier transform of the projectile density
  TAComplex *fRhoT; ///< the 2-D Fourier transform of the target density
  TAFNN *fFNN; ///< to calculate the isospin-averaged n-n scattering amplitude
  /// the resulting S-matrix
  TAComplex *fSMatrix;
  /// the fitting results
  TAComplex *fAlphaj;
  doubel *fBetaj;
};

#endif
