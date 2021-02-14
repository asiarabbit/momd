/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TAHPairing.h
  \class TAHPairing
  \brief Hamiltonian for pairing model, H = -GP+P-
  P=\sum_k{a_k^+*a_{-k}^+}, P-=(P+)+=\sum_k{a_{-k}*a_k}
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2021/02/04
  \date Last modified: 2021/02/14, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include "TAHPairing.h"
#include "TAManyBodySD.h"
#include "TAManyBodySDList.h"
#include "TAManyBodySDManager.h"
#include "TADiagonalize.h"
#include "TASPState.h"
#include "TABit.h"
#include <catch2/catch.hpp>

TAHPairing::TAHPairing(const string &configFile) : TAHamiltonian(configFile){
  fG = 1.; fDelta = 0.;
  fMatrix = nullptr;
}
TAHPairing::~TAHPairing(){}

/// \retval calculate and return the 1-body part (t+u) of the ME for H
/// here's an example for the pairing model H = -GP+P-
double TAHPairing::MatrixElement1N(int rr, int cc){
  // Slater-Gondon rules: if rr and cc have more than 1 sps different, return 0. //
  TABit br((*fMBSDListM)[rr] ^ (*fMBSDListM)[cc]);
  if(br.count() > 2) return 0.;
  // the \sum_n{delta*n*N_n} term //
  if(rr != cc) return 0.;

  int n1 = fMBSDListM->GetSPState(rr, 0)->GetN();
  int n2 = fMBSDListM->GetSPState(rr, 2)->GetN();
  return 2.*fDelta*(n1+n2);
} // end of member function MatrixElement1N

/// \retval calculate and return the 2-body part v(r1, r2) of the ME for H
/// here's an example for the pairing model H = -GP+P-
double TAHPairing::MatrixElement2N(int rr, int cc){
  // Slater-Gondon rules: if rr and cc have more than 2 sps's different, return 0. //
  TABit br((*fMBSDListM)[rr] ^ (*fMBSDListM)[cc]);
  if(br.count() > 4) return 0.;

  double me = 0., phase;
  // -GP_+P_- //
  for(int k = 0; k < fNSPState; k += 2){ // loop over pairs, so k incremented by 2
    for(int j = 0; j < fNSPState; j += 2){ // loop over pairs, so j incremented by 2
      // Integral(rr, p, q, s, r, cc): <rr|a+_p*a+_q * a_s*a_r|cc>
      // i.e. the order does not change in Integral(..), as in the operators //
      phase = fMBSDListM->Integral(rr, k, k+1, j+1, j, cc);
      if(0. == phase) continue;
      me += phase;
    } // end for over creation operators a+_q
  } // end for over creation operators a+_p
  return me * fG;
} // end member function MatrixElement2N

double TAHPairing::MatrixElement3N(int rr, int cc){
  return 0.;

  // Slater-Gondon rules: if rr and cc have more than 3 sps's different, return 0. //
  TABit br((*fMBSDListM)[rr] ^ (*fMBSDListM)[cc]);
  if(br.count() > 6) return 0.;

  return 0.;
}

/// construct the Hamiltonian matrix
matrix *TAHPairing::Matrix(){
  if(fMatrix) return fMatrix;

  const int n = GetNBasis();
  fMatrix = new matrix(n, n);
  for(int i = 0; i < n; i++) for(int j = 0; j <= i; j++){
    (*fMatrix)[j][i] = (*fMatrix)[i][j] = MatrixElement(i, j);
  }

  return fMatrix;
} // end of member function Matrix

TEST_CASE("Pairing Hamiltonian", "[pair]"){
  TAHPairing *pair = new TAHPairing("../config/input.txt");
  matrix *H = pair->Matrix();
  const int n = pair->GetNBasis(); // number of Slater determinants
  CHECK((*H)[0][0] == Approx(2.).epsilon(1e-8));
  CHECK((*H)[1][4] == Approx(0.).margin(1e-8));
  CHECK((*H)[3][5] == Approx(1.).epsilon(1e-8));
  double d[n]; matrix v(n,n);
  TADiagonalize::JacobiSort(*H, n, d, v);
  for(int i = 0; i < n; i++) std::cout << d[i] << std::endl;
  CHECK(d[0] == Approx(6.).epsilon(1e-8));
  CHECK(d[1] == Approx(2.).epsilon(1e-8));
  CHECK(d[2] == Approx(2.).epsilon(1e-8));
  TATreeCol::CloseFile();
}
