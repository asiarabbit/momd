/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAHF.cxx
  \class TAHF
  \brief Hartree-Fock algorithm to calculate single-particle states assuming
  effective mean-field imposed by all other bystanding fermions. The solved
  single-particle states can be used as input to FCI and MOMD calculations.
  The output HF 1-body and two-body matrix element would only be calculated
  upon needed, and would only be calculated once.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/11/18
  \date Last modified: 2020/12/13 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <algorithm>
#include <cstring>
#include <catch2/catch.hpp>
#include "TAException.h"
#include "TADiagonalize.h"
#include "TAMath.h"
#include "TAHF.h"

using std::swap;

const double TAHF::DUMB = -9999.;

TAHF::TAHF(int nParticle, int nSPBasis){
  if(nParticle <= 0)
    TAException::Error("TAHF", "TAHF: illgal number of particles: %d", nParticle);
  if(nSPBasis < nParticle)
    TAException::Error("TAHF", "TAHF: number of SP states less than number of particles");

  fNSPBasis = nSPBasis;
  fNParticle = nParticle;
  const int nb = nSPBasis;
  const int nv = nb*(nb-1)/2; // nv*nv: dimension of twobody matrix
  // assign memory to input one- and two- body operators
  fInputOneBody = new matrixSym(nb); fInputOneBody->Initialize(DUMB);
  fInputTwoBody = new matrixSym(nv); fInputTwoBody->Initialize(DUMB);
  fHFOneBody = new matrixSym(nb);    fHFOneBody->Initialize(DUMB);
  fHFTwoBody = new matrixSym(nv);    fHFTwoBody->Initialize(DUMB);
  // allocate memory to the HF SP vector and HF SP energy
  fHF = new matrix(nb, nb);
  fHFSPEnergy = new double[nb]{};
} // end of the construction operator
TAHF::~TAHF(){
  delete fInputOneBody;  fInputOneBody = nullptr;
  delete fInputTwoBody;  fInputTwoBody = nullptr;
  delete fHFOneBody;     fHFOneBody = nullptr;
  delete fHFTwoBody;     fHFTwoBody = nullptr;
  delete [] fHFSPEnergy; fHFSPEnergy = nullptr;
  delete fHF;            fHF = nullptr;
} // end of the destruction operator

/// assign fInputOneBody matrix, user-specific
///\retval (*fInputOneBody)[alpha][beta] = h0
void TAHF::SetInputOneBodyH(double h0, int alpha, int beta){
  if(alpha < 0 || alpha >= fNSPBasis || beta < 0 || beta >= fNSPBasis)
    TAException::Error("TAHF", "SetInputOneBodyH: array index out of border.");
  if(alpha < beta) swap(alpha, beta); // so that alpha >= beta

  if(DUMB == (*fInputOneBody)[alpha][beta]) (*fInputOneBody)[alpha][beta] = h0;
  else TAException::Warn("TAHF",
    "SetInputOneBodyH: <%d|h0|%d> already assigned.", alpha, beta);
} // end of member function SetInputOneBodyH

// transform the matrix element to standard form to accommodate symmetry imposed
// by hermiticity and invariance of V over particle-exchange
// the sign change is due to antisymmetrization of the matrix elements
void TAHF::indexSym(int alpha, int gamma, int beta, int delta, int &p, int &q, double &sign){
  sign = 1.;
  if(alpha < gamma){ swap(alpha, gamma); sign *= -1.; } // so that alpha > gamma
  if(beta < delta){ swap(beta, delta); sign *= -1.; }// so that beta > delta
  p = alpha*(alpha-1)/2 + gamma; q = beta*(beta-1)/2 + delta;
  if(p < q) swap(p, q); // so that p >= q, no sign change due to hermiticity of V
}
// note that the matrix element here is antisymmetrized, i.e.,
// = <alpha gamma|V|beta delta>_{AS}= v - vm
// = <alpha gamma|V|beta delta> - <alpha gamma|V|delta beta>
// v and vm
void TAHF::SetInputTwoBodyH(double v, double vm, int alpha, int gamma,
    int beta, int delta){
  if(alpha < 0 || alpha >= fNSPBasis || beta < 0 || beta >= fNSPBasis ||
     beta < 0 || beta >= fNSPBasis || delta < 0 || delta >= fNSPBasis)
    TAException::Error("TAHF", "SetInputTwoBodyH: array index out of border.");
  if(alpha == gamma || beta == delta){
    TAException::Warn("TAHF", "SetInputTwoBodyH: two particle at the same state!");
    return; // Pauli's exclusion principle
  }
  // <alpha,gamma|v|beta,delta> -> <p|v|q>
  // using symmetry and antisymmetrization of v to save space
  int p, q; double sign; indexSym(alpha, gamma, beta, delta, p, q, sign);

  // get the subscript finally //
  if(DUMB == (*fInputTwoBody)[p][q]) (*fInputTwoBody)[p][q] = (v - vm) * sign;
  else TAException::Warn("TAHF",
    "SetInputTwoBodyH: designated matrix element already assigned.\
alpha: %d, gamma: %d, gamma: %d, delta: %d", alpha, gamma, beta, delta);
} // end of member function SetInputTwoBodyH
// retrieve the one- and two-body part of the hamiltonian in original basis
// with symmetry considered
double TAHF::OneBodyH(int alpha, int beta){
  if(alpha < 0 || alpha >= fNSPBasis || beta < 0 || beta >= fNSPBasis)
    TAException::Error("TAHF", "OneBodyH: array index out of border.");

  if(alpha < beta) swap(alpha, beta);
  if(DUMB == (*fInputOneBody)[alpha][beta]) TAException::Warn("TAHF",
    "OneBodyH: Required element not assignd yet. alpha: %d, beta: %d", alpha, beta);
  return (*fInputOneBody)[alpha][beta];
} // end of member function OneBodyH

double TAHF::TwoBodyH(int alpha, int gamma, int beta, int delta){
  if(alpha < 0 || alpha >= fNSPBasis || beta < 0 || beta >= fNSPBasis ||
     beta < 0 || beta >= fNSPBasis || delta < 0 || delta >= fNSPBasis)
    TAException::Error("TAHF", "TwoBodyH: array index out of border.");
  if(alpha == gamma || beta == delta){
    TAException::Warn("TAHF", "TwoBodyH: two particle at the same state!");
    return 0.; // Pauli's exclusion principle
  }
  // <alpha,gamma|v|beta,delta> -> <p|v|q>
  // using symmetry and antisymmetrization of v to save space
  int p, q; double sign; indexSym(alpha, gamma, beta, delta, p, q, sign);
  if(DUMB == (*fInputOneBody)[p][q])
    TAException::Warn("TAHF", "TwoBodyH: designated matrix element not assigned yet.\
      alpha: %d, gamma: %d, gamma: %d, delta: %d", alpha, gamma, beta, delta);
  return (*fInputOneBody)[p][q] * sign;
} // end of member function TwoBodyH


/// the Hartree-Fock iteration
void TAHF::HF(){
  // so that HF() is only called once //
  static bool isCalled = false;
  if(isCalled) return;
  isCalled = true;

  static const int NITR = 100; // maximum number of iterations
  static const double EPS = 1.e-5; // precision in HF s.p. energy, unit: MeV
  const int np = GetNParticle(), nb = GetNSPBasis();
  Initialize(); // initialize fHF
  matrix rho(nb, nb); // the density matrix to build HF effective mean field
  double ehf[nb]{}; // the Hartree-Fock energy of the last iteration
  // here the iteration begins //
  for(int iter = 0; iter < NITR; iter++){
    // 1. calculate the density matrix //
    rho = 0.; // initialize the density matrix
    for(int gamma = nb; gamma--;) for(int delta = nb; delta--;) for(int j = np; j--;)
      rho[gamma][delta] += (*fHF)[gamma][j]*(*fHF)[delta][j]; // C_gammaj*C^*_deltaj
    // 2. construct the HF matrix (Fock matrix), note that fHF is a symmetric matrix //
    for(int alpha = nb; alpha--;) for(int beta = alpha+1; beta--;){ // the lower triangle
      (*fHF)[alpha][beta] = OneBodyH(alpha, beta);
      // build the effective mean field //
      for(int gamma = nb; gamma--;) for(int delta = nb; delta--;)
        (*fHF)[alpha][beta] += rho[delta][gamma] * TwoBodyH(alpha, gamma, beta, delta);
      if(alpha != beta) (*fHF)[beta][alpha] = (*fHF)[alpha][beta];
    } // end of for over alpha and beta
    // 3. diagonalize fHF //
    memcpy(ehf, fHFSPEnergy, sizeof(ehf));
    // fHF replaced by its own eigenvectors here
    TADiagonalize::TridiagQLSort(*fHF, nb, fHFSPEnergy);
    // 4. test convergence //
    for(int i = nb; i--;) ehf[i] -= fHFSPEnergy[i];
    if(TAMath::infNorm(nb, ehf) < EPS) return;
  } // end of the main iterations
  double eps = TAMath::infNorm(nb, ehf);
  TAException::Warn("TAHF", "HF: Too many iterations occurred. iter: %d, eps: %f", NITR, eps);
  return; // never gets here
} // end of member function HF

// retrieve the one- and two-body part of the hamiltonian in the calculated HF basis
double TAHF::OneBodyHF(int p, int q){
  if(p < 0 || p >= fNSPBasis || q < 0 || q >= fNSPBasis)
    TAException::Error("TAHF", "OneBodyHF: array index out of border.");
  if(p < q) swap(p, q);

  if(DUMB == (*fHFOneBody)[p][q]){
    (*fInputOneBody)[p][q] = 0.; HF(); // HF() has to be called at least once
    for(int alpha = fNSPBasis; alpha--;) for(int beta = fNSPBasis; beta--;)
      (*fInputOneBody)[p][q] += (*fHF)[alpha][p]*OneBodyH(alpha, beta)*(*fHF)[beta][q];
  } // end if
  return (*fHFOneBody)[p][q];
} // end of member function OneBodyHF
double TAHF::TwoBodyHF(int p, int q, int r, int s){
  if(p < 0 || p >= fNSPBasis || q < 0 || q >= fNSPBasis ||
     r < 0 || r >= fNSPBasis || s < 0 || s >= fNSPBasis)
    TAException::Error("TAHF", "TwoBodyHF: array index out of border.");
  if(p == q || r == s){
    TAException::Warn("TAHF", "TwoBodyHF: two particle at the same state!");
    return 0.; // Pauli's exclusion principle
  }
  // <alpha,gamma|v|beta,delta> -> <p|v|q>
  // using symmetry and antisymmetrization of v to save space
  int pp, qq; double sign; indexSym(p, q, r, s, pp, qq, sign);
  if(DUMB == (*fHFTwoBody)[pp][qq]){
    (*fHFTwoBody)[pp][qq] = 0.; HF(); // HF() has to be called at least once
    for(int alpha = fNSPBasis; alpha--;) for(int beta  = fNSPBasis;  beta--;)
    for(int gamma = fNSPBasis; gamma--;) for(int delta = fNSPBasis; delta--;)
      (*fHFTwoBody)[pp][qq] += (*fHF)[alpha][p]*(*fHF)[gamma][q]*
        TwoBodyH(alpha, gamma, beta, delta)*(*fHF)[beta][r]*(*fHF)[delta][s];
  } // end if
  return (*fHFTwoBody)[pp][qq] * sign;
} // end of member function TwoBodyHF

// initialize fHF to initiate the iteration
void TAHF::Initialize(){
  *fHF = 1.; // initilize to delta_{alpha,i}
} // end of member function Initialize
