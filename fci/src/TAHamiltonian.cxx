/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TAHamiltonian.cxx
  \class TAHamiltonian
  \brief The hamiltonian itself, period. It is responsible for storing the
  formula (int string form), the operator (creation and annhilation) form, and
  the matrix form. When in matrix form, the basis of the representation would
  also be specified as a member of this class. We deem it appropriate to make
  this a singleton class.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/09
  \date Last modified: 2020/02/28, by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include "TAManyBodySD.h"
#include "TAHamiltonian.h"
#include "TAManyBodySDList.h"
#include "TAManyBodySDManager.h"
#include "TASingleParticleState.h"
#include "TASingleParticleStateManager.h"


TAHamiltonian *TAHamiltonian::kInstance = nullptr;

TAHamiltonian::TAHamiltonian() : fCoe1N(0), fCoe2N(0), fCoe3N(0),
  fMBSDListM(0), fMatrix(0), fNSPState(0), fNMBSD(0){
  // prepare the basis of the representation //
  TAManyBodySDManager *mbsdManager = TAManyBodySDManager::Instance();
  mbsdManager->MSchemeGo(); // generate many-body basis
  SetMBSDListM(mbsdManager->GetMBSDListM());
  fNMBSD = fMBSDListM->GetNBasis();
  fNSPState = TASingleParticleStateManager::Instance()->GetNSPState();
} // end of the constructor

TAHamiltonian::~TAHamiltonian(){
  if(fCoe1N){ delete fCoe1N; fCoe1N = nullptr; }
  if(fCoe2N){ delete fCoe2N; fCoe2N = nullptr; }
  if(fCoe3N){ delete fCoe3N; fCoe3N = nullptr; }
  if(fMatrix){ delete fMatrix; fMatrix = nullptr; }
} // end of the destructor

TAHamiltonian *TAHamiltonian::Instance(){
  if(!kInstance) kInstance = new TAHamiltonian();
  return kInstance;
}

/// \retval calculate and return the matrix form of the hamiltonian
matrix &TAHamiltonian::Matrix(){
  if(fMatrix && !fMatrix->IsEmpty()){
    return *fMatrix;
  }
  if(!fCoe1N){
    TAException::Error("TAHamiltonian",
      "Matrix: 1-body operator coefficient matrix not assigned.");
  }
  if(!fCoe2N){
    TAException::Error("TAHamiltonian",
      "Matrix: 2-body operator coefficient matrix not assigned.");
  }
  if(!fMBSDListM){
    TAException::Error("TAHamiltonian",
      "Matrix: The many-body basis not assigned.");
  }
  if(!fNSPState){
    TAException::Error("TAHamiltonian", "Matrix: fNSPState not assigned.");
  }
  if(!fNMBSD){
    TAException::Error("TAHamiltonian", "Matrix: fNMBSD is 0. Not assigned?");
  }

  // loop to generate each matrix element for the hamiltonian //
  if(fMatrix){ delete fMatrix; fMatrix = nullptr; }
  fMatrix = new matrix(fNMBSD, fNMBSD); // allot memery to a nxn matrix
  // initialize to a specific initial value //
  for(int i = fNMBSD; i--;) for(int j = fNMBSD; j--;) (*fMatrix)[i][j] = -9999.;

  for(int rr = 0; rr < fNMBSD; rr++){
    for(int cc = 0; cc < fNMBSD; cc++){
      MatrixElement(rr, cc); // assign matrix element H[i][j]
    } // end for over columns
  } // end for over rows

  return *fMatrix;
} // end of the member function Matrix

/// assign the matrix element (*fMatrix)[i][j]
/// \param r: row, c: column
void TAHamiltonian::MatrixElement(int rr, int cc){
  // so that the matrix is symmetric //
  if((*fMatrix)[rr][cc] != -9999.) return;
  if((*fMatrix)[cc][rr] != -9999.){
    (*fMatrix)[rr][cc] = (*fMatrix)[cc][rr];
    return;
  }
  (*fMatrix)[rr][cc] = MatrixElement1N(rr, cc) +
    MatrixElement2N(rr, cc) + MatrixElement3N(rr, cc);
} // end of member function MatrixElement

/// \retval calculate and return the 1-body part (t+u) of the ME for H
double TAHamiltonian::MatrixElement1N(int rr, int cc){
  double me = 0., phase, force;
  for(int p = 0; p < fNSPState; p++){
    for(int q = 0; q < fNSPState; q++){
      /// Integral(rr, p, q, cc): <rr|a+_p * a_q|cc>
      if(!(phase = fMBSDListM->Integral(rr, p, q, cc)) ||
         !(force = (*fCoe1N)[p][q]) ) continue;
      me += force * phase;
    } // end for over annhilation operators a_q
  } // end for over creation operators a+_p
  return me;
} // end member function MatrixElement1N

/// \retval calculate and return the 2-body part v(r1, r2) of the ME for H
double TAHamiltonian::MatrixElement2N(int rr, int cc){
  double me = 0., phase, force;
  if(!fCoe2N) return 0.; // 2N force is not assigned
  for(int p = 0; p < fNSPState; p++){
    for(int q = 0; q < fNSPState; q++){
      if(q == p) continue; // Pauli's exclusion principle
      for(int r = 0; r < fNSPState; r++){ // -------------------------------
        for(int s = 0; s < fNSPState; s++){
          if(s == r) continue; // Pauli's exclusion principle
          /// Integral(rr, p, q, r, s, cc): <rr|a+_p*a+_q * a_s*a_r|cc>
          if(!(phase = fMBSDListM->Integral(rr, p, q, s, r, cc)) ||
             !(force = (*fCoe2N)[p][q][r][s]) ) continue;
//          cout << "<rr|p+ q+ s r|cc>" << endl; // DEBUG
//          printf("<%d|%d+ %d+ %d %d|%d>\n", rr, p, q, s, r, cc); // DEBUG
//          cout << "force: " << force << " phase: " << phase << endl; // DEBUG
//          getchar(); // DEBUG
          me += force * phase;
        } // end for over annhilation operators a_s
      } // end for over annhilation operators a_r // -----------------------
    } // end for over creation operators a+_q
  } // end for over creation operators a+_p
//  TAException::Debug("Mark1", "rr: %d, cc: %d, me: %f", rr, cc, me); // DEBUG
  return me / 4.; // 4 = (2!)^2
} // end member function MatrixElement2N

/// \retval calculate and return the 3-body part v(r1,r2,r3) of the ME for H
double TAHamiltonian::MatrixElement3N(int rr, int cc){
  if(!fCoe3N) return 0.; // 3N force is not needed
  double me = 0., phase, force;
  for(int p = 0; p < fNSPState; p++){
    for(int q = 0; q < fNSPState; q++){
      if(p == q) continue; // Pauli's exclusion principle
      for(int r = 0; r < fNSPState; r++){
        if(r == p || r == q) continue; // Pauli's exclusion principle
        for(int s = 0; s < fNSPState; s++){ // -----------------------------
          for(int t = 0; t < fNSPState; t++){
            if(t == s) continue; // Pauli's exclusion principle
            for(int u = 0; u < fNSPState; u++){
              if(u == s || u == t) continue; // Pauli's exclusion principle
              /// Integral(rr, p, q, r, s, t, u, cc):
              /// <rr|a+_p*a+_q*a+_r * a_u*a_t*a_s|cc>
              if(!(phase = fMBSDListM->Integral(rr, p, q, r, u, t, s, cc)) ||
                 !(force = (*fCoe3N)[p][q][r][s][t][u]) ) continue;
              me += force * phase;
            } // end for over annhilation operators a_u
          } // end for over annhilation operators a_t
        } // end for over annhilation operators a_s // ---------------------
      } // end for over creation operators a+_r
    } // end for over creation operators a+_q
  } // end for over creation operators a+_p
  return me / 36.; // 36 = (3!)^2
} // end member function MatrixElement3N

void TAHamiltonian::SetCoe1N(const matrix &coe1N){
  if(fCoe1N){ delete fCoe1N; fCoe1N = nullptr; }
  fCoe1N = new matrix(coe1N);
} // end of member function SetCoe1N
void TAHamiltonian::SetCoe2N(const TAMatrix4D &coe2N){
  if(fCoe2N){ delete fCoe2N; fCoe2N = nullptr; }
  fCoe2N = new TAMatrix4D(coe2N);
} // end of member function SetCoe2N
void TAHamiltonian::SetCoe3N(const TAMatrix6D &coe3N){
  if(fCoe3N){ delete fCoe3N; fCoe3N = nullptr; }
  fCoe3N = new TAMatrix6D(coe3N);
} // end of member function SetCoe3N
void TAHamiltonian::SetMBSDListM(TAManyBodySDList *mbsd){
  if(!mbsd){
    TAException::Error("TAHamiltonian", "SetMBSDListM: Input pointer is null.");
  }
  fMBSDListM = mbsd;
} // end of member function SetMBSDListM

// so that this class could undergo a debugging test //
// the coefficients of 1, 2 and 3-body force are all set to 1 //
void TAHamiltonian::InitializeCoefficient(){
  if(!fNSPState){
    TAException::Error("TAHamiltonian",
      "InitializeCoefficient: fNSPState not set yet.");
  }

  if(fCoe1N){ delete fCoe1N; fCoe1N = nullptr; }
  if(fCoe2N){ delete fCoe2N; fCoe2N = nullptr; }
  if(fCoe3N){ delete fCoe3N; fCoe3N = nullptr; }

  // initialize fCoe1N //
  fCoe1N = new matrix(fNSPState, fNSPState);
  // initialze fCoe1
  const vector<TASingleParticleState *> &spsVec =
    TASingleParticleStateManager::Instance()->GetSPStateVec();
  for(int i = 0; i < fNSPState; i++){
    for(int j = 0; j < fNSPState; j++){
      if(i != j) (*fCoe1N)[i][j] = 0.;
      else (*fCoe1N)[i][j] = spsVec[i]->GetEnergy();
    } // end for over columns
  } // end for over rows

  // initialize nucleon-correlation force coefficients //
  matrix tmp(fNSPState, fNSPState); // all elements would be set to 1 //
  tmp = 0.;
  for(int i = 0; i < fNSPState; i++){
    for(int j = i + 1; j < fNSPState; j++){
      // apply antisymmetrization //
      tmp[i][j] = 1.;
      tmp[j][i] = -1.;
    } // end for over j
  } // end for over i
  // initialize fCoe2N //
  matrix tmpMinus(tmp); tmp *= -1.;
  fCoe2N = new TAMatrix4D(fNSPState, fNSPState);
  for(int i = 0; i < fNSPState; i++){
    for(int j = i; j < fNSPState; j++){
      // apply antisymmetrization //
      if(i == j){
        (*fCoe2N)[i][j].Resize(fNSPState, fNSPState);
        (*fCoe2N)[i][j] = 0.;
      }
      else{
        (*fCoe2N)[i][j] = tmp;
        (*fCoe2N)[j][i] = tmpMinus;
      }
    } // end for over j
  } // end for over i
  return; // DEBUG
  // initialize fCoe3N //
  fCoe3N = new TAMatrix6D(fNSPState, fNSPState);
  for(int i = 0; i < fNSPState; i++){
    for(int j = 0; j < fNSPState; j++){
      // apply space first
      (*fCoe3N)[i][j] = (*fCoe2N);
    } // end for over j
  } // end for over i
  // then implement the antisymmetrization //
  TABit bitL, bitR; // bit representation of a many-body Slater determinant
  for(int p = 0; p < fNSPState; p++){
    for(int q = 0; q < fNSPState; q++){
      for(int r = 0; r < fNSPState; r++){
        for(int s = 0; s < fNSPState; s++){
          for(int t = 0; t < fNSPState; t++){
            for(int u = 0; u < fNSPState; u++){
              /// <bitL|a+_p*a+_q*a+_r * a_s*a_t*a_u|bitR>
              bitL.Reset(); bitR.Reset();
              (*fCoe3N)[p][q][r][s][t][u] = 1. *
                bitL.Create(r).Create(q).Create(p).GetPhase() *
                bitR.Create(u).Create(t).Create(s).GetPhase();
            } // end for over annhilation operators a_u
          } // end for over annhilation operators a_t
        } // end for over annhilation operators a_s
      } // end for over creation operators a+_r
    } // end for over creation operators a+_q
  } // end for over creation operators a+_p
} // end of member function InitializeCoefficient
