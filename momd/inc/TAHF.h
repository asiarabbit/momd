/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAHF.h
  \class TAHF
  \brief Hartree-Fock algorithm to calculate single-particle states assuming
  effective mean-field imposed by all other bystanding fermions. The solved
  single-particle states can be used as input to FCI and MOMD calculations.
  Note that user has to select original s.p basis and one-body and two-body
  matrix elements of the Hamiltonian, which also is upon the choice of the
  user.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/11/17
  \date Last modified: 2020/11/17 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/


#ifndef _TAHF_h_
#define _TAHF_h_

#include "TAMatrix.h"
#include "TASymMatrix.h"

typedef TASymMatrix<double> matrixSym;

class TAHF{
public:
  TAHF(int nParticle, int nSPBasis);
  virtual ~TAHF();

  int GetNParticle() const{ return fNParticle; }
  int GetNSPBasis() const{ return fNSPBasis; }
  /// assign fInputOneBody matrix, user-specific
  ///\retval fInputOneBody[alpha][beta] = h0
  virtual void SetInputOneBodyH(double h0, int alpha, int beta);

  // transform the matrix element to standard form to accommodate symmetry imposed
  // by hermiticity and invariance of V over particle-exchange
  // the sign change is due to antisymmetrization of the matrix elements
  static void indexSym(int alpha, int gamma, int beta, int delta, int &p,
    int &q, double &sign);
  /// assign fInputTwoBody matrix, user-specific
  ///\retval fInputTwoBody[alpha][gamma][beta][delta] = v
  virtual void SetInputTwoBodyH(double v, double vm, int alpha, int gamma,
    int beta, int delta);
  // retrieve the one- and two-body part of the hamiltonian in original basis
  virtual double OneBodyH(int alpha, int beta);
  virtual double TwoBodyH(int alpha, int gamma, int beta, int delta);
  virtual void HF(); ///< the Hartree-Fock iteration
  virtual void Initialize(); ///< initialize fHFSPVector to initiate the iteration
  // retrieve the one- and two-body part of the hamiltonian in the calculated HF basis
  double OneBodyHF(int p, int q);
  double TwoBodyHF(int p, int q, int r, int s);

  // marking value for those that are not assigned
  static const double DUMB;

private:
  /// THE INPUT ///
  int fNParticle;
  int fNSPBasis;
  matrixSym *fInputOneBody; ///< the one-body part of H in original SP basis
  matrixSym *fInputTwoBody; ///< the two-body part of H in original SP basis

  /// THE OUTPUT ///
  /// note that upon output, the one- and two-body part of H should undergo a
  /// unitary transform since we'll update the representation to HF basis in FCI
  /// calculations
  matrixSym *fHFOneBody; ///< the one-body part of H in HF SP basis
  matrixSym *fHFTwoBody; ///< the two-body part of H in HF SP basis
  /// the eigenvector and eigenvalue of matrix fHF
  matrix *fHF; ///< the Fock matrix, together it stores its own eigenvectors
  double *fHFSPEnergy; ///< the solved HF single-particle energy, i.e. fHF's eigenvalues
};

#endif
