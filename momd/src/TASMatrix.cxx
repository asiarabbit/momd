/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASMatrix.cxx
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


TASMatrix::TASMatrix(){}
virtual ~TASMatrix::TASMatrix(){}

/// calculate and output the S-matrix
void TASMatrix::SMatrix(TAComplex *smat){}
/// fit the S-matrix to an expansion of Gaussians
/// the results are stored in member arrays alphaj and betaj
void TASMatrix::GaussianFit(double *alphaj, double *betaj){}

void TASMatrix::GetAlphaj(double *alphaj){}
void TASMatrix::GetBetaj(double *betaj){}
