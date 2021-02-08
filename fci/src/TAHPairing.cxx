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
#include "TAManyBodySD.h"
#include "TAManyBodySDList.h"
#include "TASingleParticleState.h"

TAHPairing::TAHPairing(const string &configFile) : TAHamiltonian(configFile){
  fG = 1.; fDelta = 0.;
  fMBSDListM->Pairing();
}
TAHPairing::~TAHPairing(){}


/// \retval calculate and return the 1-body part (t+u) of the ME for H
/// here's an example for the pairing model H = -GP+P-
double TAHPairing::MatrixElement1N(int rr, int cc){
  // Slater-Gondon rules: if rr and cc have more than 1 sps different, return 0. //
  TAManyBodySD *sdr = (*fMBSDListM)[rr], *sdc = (*fMBSDListM)[cc];
  auto br = sdr->Bit(); br ^= sdc->Bit();
  if(br.count() > 2) return 0.;

  // the \sum_n{delta*n*N_n} term //
  if(rr != cc) return 0.;
  int n1 = (*sdr)[0]->GetN(), n2 = (*sdc)[2]->GetN();
  return 2.*fDelta*(n1+n2);
} // end of member function MatrixElement1N

/// \retval calculate and return the 2-body part v(r1, r2) of the ME for H
/// here's an example for the pairing model H = -GP+P-
double TAHPairing::MatrixElement2N(int rr, int cc){
  // Slater-Gondon rules: if rr and cc have more than 2 sps's different, return 0. //
  TAManyBodySD *sdr = (*fMBSDListM)[rr], *sdc = (*fMBSDListM)[cc];
  TABit br(sdr->Bit()); br ^= sdc->Bit();
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
  TAManyBodySD *sdr = (*fMBSDListM)[rr], *sdc = (*fMBSDListM)[cc];
  auto br = sdr->Bit(); br ^= sdc->Bit();
  if(br.count() > 6) return 0.;

  return 0.;
}
