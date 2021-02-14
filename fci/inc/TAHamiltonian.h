/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TAHamiltonian.h
  \class TAHamiltonian
  \brief The Hamiltonian matrix class, to generate and manage the many-body Slater
  determinants, and also interfaces the matrix element calculation methods. Note
  that this is an abstract base class. Calculation of the matrix elements should
  be specifically implemented in its sub-class.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/09
  \date Last modified: 2021/02/04, by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAHamiltonian_h_
#define _TAHamiltonian_h_

#include <string>
#include "TAMatrix.h"
#include "TASparseVec.h"

using std::string;
class TAManyBodySDList;
class TAHF; // to provide Hartree-Fock basis

class TAHamiltonian{
public:
  /// the constructor, made private to be identified as a singleton class
  TAHamiltonian(const string &configFile);
  virtual ~TAHamiltonian();
  /// \retval calculate and return the matrix form of the hamiltonian
  unsigned long GetNBasis() const;

  /// assign the matrix element (*fMatrix)[i][j]
  virtual double MatrixElement(int row, int column);
  /// the following are supposed to be user-specific, i.e. the definitions of the
  /// following methods are Hamiltonian-dependent, and must defined by users.
  virtual double MatrixElement1N(int rr, int cc) = 0;
  virtual double MatrixElement2N(int rr, int cc) = 0;
  virtual double MatrixElement3N(int rr, int cc) = 0;
  void DotProduct(TASparseVec &v, TASparseVec &r); /// r=H.v

protected:
  /// M-scheme many-body SD list, to define the representation
  TAManyBodySDList *fMBSDListM; ///< \NOTE its memory doesn't need to be freed
  int fNSPState; ///< number of single particle states
  int fNParticle; ///< number of particles
};

#endif
