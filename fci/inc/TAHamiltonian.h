/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TAHamiltonian.h
  \class TAHamiltonian
  \brief The hamiltonian itself, period. It is responsible for storing the
  formula (int string form), the operator (creation and annhilation) form, and
  the matrix form. When in matrix form, the basis of the representation would
  also be specified as a member of this class. We deem it appropriate to make
  this a singleton class.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/09
  \date Last modified: 2020/02/12, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAHamiltonian_h_
#define _TAHamiltonian_h_

#include <string>
#include "TAMatrix.h"

typedef TAMatrix<matrix> TAMatrix4D; /// 4-th order tensor
typedef TAMatrix<TAMatrix4D> TAMatrix6D; /// 6-th order tensor

class TAManyBodySDList;

using std::string;

class TAHamiltonian{
public:
  virtual ~TAHamiltonian();
  static TAHamiltonian *Instance();
  /// \retval return the specific formula of the hamiltonian
  const char *Formula() const{ return fFormula.c_str(); }
  /// \retval calculate and return the matrix form of the hamiltonian
  matrix &Matrix();
  vec_t<double> &operator[](int i){ return Matrix()[i]; }

  void SetCoe1N(const matrix &coe1N);
  void SetCoe2N(const TAMatrix4D &coe2N);
  void SetCoe3N(const TAMatrix6D &coe3N);
  void SetMBSDListM(TAManyBodySDList *mbsd);

  /// the coefficients of 1, 2 and 3-body force are all set to 1 //
  /// so that this class could undergo a debugging test //
  void InitializeCoefficient();

private:
  /// the constructor, made private to be identified as a singleton class
  TAHamiltonian();
  /// assign the matrix element (*fMatrix)[i][j]
  void MatrixElement(int row, int column);
  /// \retval calculate and return the 1-body part (t+u) of the ME for H
  double MatrixElement1N(int row, int column);
  /// \retval calculate and return the 2-body part v(r1, r2) of the ME for H
  double MatrixElement2N(int row, int column);
  /// \retval calculate and return the 3-body part v(r1,r2,r3) of the ME for H
  double MatrixElement3N(int row, int column);

  static TAHamiltonian *kInstance;
  /// \NOTE note that all these coeffiicients are supoosed to be user input
  matrix *fCoe1N; ///< coefficients for the 1-N part of H: <p|t+u|q>
  TAMatrix4D *fCoe2N; ///< coefficients for the 2-N part of H: <pq|v|rs>
  TAMatrix6D *fCoe3N; ///< coefficients for the 3-N part of H: <pqr|v_3|stu>
  /// M-scheme many-body SD list, to define the representation
  TAManyBodySDList *fMBSDListM; ///< \NOTE its memory doesn't need to be freed
  matrix *fMatrix; ///< the hamiltonian matrix in fMBSDListM basis
  int fNSPState; ///< number of single particle states
  int fNMBSD; ///< number of many-body Slater determinants in fMBSDListM
  string fFormula;
};

#endif
