/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TALanczos.h
  \class TALanczos
  \brief Lanczos algorithm for solving large sparse matrix iteratively. Transcribbed
  from momd/TADiagonalize::Lanczos and momd/TADiagonalize::LanczosPurge, but with
  new classes representing large sparse matrix using ROOT - TTree object storing
  large sparse matrix and vectors, so as to save memory. This is a static class.
  Dedicated for diagonalization problems arising from large-scale shell model calculations.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/12, Chinese Lunar New Year
  \date Last modified: 2021/02/12 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TALanczos_h_
#define _TALanczos_h_

#include "TAMatrix.h"

class TAHamiltonian;
class TASparseMatrix;

class TALanczos{
public:
  TALanczos(){}
  virtual ~TALanczos(){}

  /// the corresponding eigenvectors of huge sparse matrix A
  /// Lanczos algorithm is devised for a few largest (in modulus) eigenvalues and
  /// It projects A to a space with much lower dimension (Krylov subspace) as a tridiagonal
  /// matrix, whose eigenvalues and eigenvectors are good approximation of A's
  /// Given input matrix a[0..n-1][0..n-1], the routine returns first fiew (usually <=10)
  /// eigenvalues in d and the corresponding eigenvectors in z
  /// d should be at least of length 50, for it is drafted in the program
  static void Lanczos(TAHamiltonian &H, double *d, TASparseMatrix &z);
  /// Lanczos method featuring Ritz-pairs purging. Unwanted Ritz-pairs are deflated,
  /// leaving a smaller Krylov space, where the Lanczos orthogonalization resumes
  /// so that the unwanted Ritz-pairs are purged from the final constructed Krylov
  /// space. This would speed the convergence to the wanted Ritz-pairs compared with
  /// plain Lanczos algorithm. NOTE that we assume the largest eigenvalues are wanted
  /// Ref.: http://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf, p212
  /// Ref.: K. Wu and H. D. Simon,SIAM J. Matrix Anal. Appl., 22 (2000), pp. 602–616
  static void LanczosPurge(TAHamiltonian &H, double *d, TASparseMatrix &z);

};

#endif
