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

#include "TAHPairing.h"

TAHPairing::TAHPairing(const string &configFile) : TAHamiltonian(configFile){
  fG = 1.;
}
TAHPairing::~TAHPairing(){}


/// \retval calculate and return the 1-body part (t+u) of the ME for H
/// here's an example for the pairing model H = -GP+P-
double TAHPairing::MatrixElement1N(int rr, int cc){ return 0.; }
/// \retval calculate and return the 2-body part v(r1, r2) of the ME for H
/// here's an example for the pairing model H = -GP+P-
double TAHPairing::MatrixElement2N(int rr, int cc){
  double me = 0., phase;
  for(int k = 0; k < fNSPState; k+=2){ // loop over pairs, so k incremented by 2
    for(int j = 0; j < fNSPState; j+=2){ // loop over pairs, so j incremented by 2
      /// Integral(rr, p, q, r, s, cc): <rr|a+_p*a+_q * a_s*a_r|cc>
      phase = fMBSDListM->Integral(rr, k, k+1, j+1, j, cc);
      if(!phase) continue;
      me += phase;
    } // end for over creation operators a+_q
  } // end for over creation operators a+_p
  return me*fG; // 4 = (2!)^2
} // end member function MatrixElement2N

double TAHPairing::MatrixElement3N(int rr, int cc){
  return 0.;
}
