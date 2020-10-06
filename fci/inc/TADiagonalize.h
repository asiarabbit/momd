/**
  \file TADiagonalize.h
  \class TADiagonalize
  \brief A collection of matrix diagonalization methods.
  \author SUN Yazhou
  \date Created: 2020/09/27
  \date Last modified: 2020/09/29 by SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TADiagonalize_h_
#define _TADiagonalize_h_

#include "TAMatrix.h"

class TADiagonalize{
public:
  TADiagonalize();
  virtual ~TADiagonalize();

  /// computes all eigenvalues and eigenvectors of a real symmetric matrix
  /// a[0..n-1][0..n-1]. On output, elements of a are destroyed. d[0..n-1] returns
  /// the eigenvalues of a. v is a nxn matrix whose columns contain, on output, the
  /// normalized eigenvectors of a. The function returns the number of Jacobi rotations
  /// that were required.
  /// This implementation is transcribbed from Ref. Numerical Recipes in C, p467.
  /// for machines where underflow is not set to zero, the program has to be modified.
  static int Jacobi(matrix &a, int n, double *d, matrix &v);
  /// note that eigenvalues are not ordered in Jacobi. This routine below sorts
  /// the input eigenvalues into descending order, and rearranges the eigenvectors
  /// correspondingly. The method is straight insertion.
  static void EigenSort(double *d, matrix &v, int n);
  /// Householder reduction of a real symmetric matrix to tridiagonal form
  /// On output, a is replaced by the orthogonal matrix Q effecting the transformation
  /// d[0..n-1] returns the diagonal and e[0..n-1] the sub-diagonal, with e[0]=0
  /// This implementation is transcribbed from Numerical Recipes in C, p474
  static void Tridiagonalize(matrix &a, int n, double *d, double *e);
  /// QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors
  /// of a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix previously
  /// reduced by Tridiagonalize(...) above. On input, d[0..n-1] contains the diagonal
  /// elements of the tridiagonal matrix. On output, it returns the eigenvalues. The vector
  /// e[0..n-1] inputs the subdiagonal elements of the tridiagonal matrix, with e[0] arbitray.
  /// On output e is destroyed. When finding only the eigenvalues, several lines may be omitted,
  /// as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired,
  /// the matrix z[0..n-1][0..n-1] is input a sthe identity matrix. if the eigenvectors of
  /// a matrix that has been reduced by Tridiagonalize are required, then z is input as
  /// the matrix output by Tridiagonalize. In either case, the kth column of z returns the
  /// normalized eigenvector corresponding to d[k].
  static void TridiagQLImplicit(double *d, double *e, int n, matrix &z);
  /// An initial attempt for implementation of Lanczos algorithm ///
  /// Lanczos algorithm is devised for a few largest (in modulus) eigenvalues and
  /// the corresponding eigenvectors of huge sparse matrix A
  /// It projects A to a space with much lower dimension (Krylov subspace) as a tridiagonal
  /// matrix, whose eigenvalues and eigenvectors are good approximation of A's
  /// Given input matrix a[0..n-1][0..n-1], the routine returns first fiew (usually <=10)
  /// eigenvalues in d and the corresponding eigenvectors in z
  static void Lanczos(matrix &a, int n, double *d, matrix &z, double *x = nullptr);
};

#endif
