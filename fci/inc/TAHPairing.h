/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TAHPairing.h
  \class TAHPairing
  \brief Hamiltonian for pairing model, H = -GP+P-
  P=\sum_k{a_k^+*a_{-k}^+}, P-=(P+)+=\sum_k{a_{-k}*a_k}
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2021/02/04
  \date Last modified: 2021/02/04, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAHPairing_h_
#define _TAHPairing_h_

#include "TAHamiltonian.h"
#include "TAManyBodySDList.h"
#include "TAMatrix.h"

class TAHPairing : public TAHamiltonian{
public:
  TAHPairing(const string &configFile);
  virtual ~TAHPairing();

  /// the following are supposed to be user-specific, i.e. the definitions of the
  /// following methods are Hamiltonian-dependent, and must defined by users.
  virtual double MatrixElement1N(int rr, int cc) override;
  virtual double MatrixElement2N(int rr, int cc) override;
  virtual double MatrixElement3N(int rr, int cc) override;
  matrix *Matrix(); ///< construct the Hamiltonian matrix

protected:
  /// H=\sum_n{delta*n*N_n} - GP_+P_-
  double fG, fDelta; ///< the interaction strength constant
  matrix *fMatrix; ///< the pairing hamiltonian matrix
};

#endif
