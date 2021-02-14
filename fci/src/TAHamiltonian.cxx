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

#include "TAManyBodySD.h"
#include "TAHamiltonian.h"
#include "TAManyBodySDList.h"
#include "TAManyBodySDManager.h"
#include "TASPState.h"
#include "TASPStateManager.h"

static const double DUMB = -9999.;

TAHamiltonian::TAHamiltonian(const string &configFile)
    : fMBSDListM(0), fNSPState(0){
  // prepare the single-particle basis set //
  // prepare the basis of the representation //
  TAManyBodySDManager *mbsdManager = TAManyBodySDManager::Instance();
  mbsdManager->LoadConfigFile(configFile);
  mbsdManager->MSchemeGo(); // generate many-body basis
  fMBSDListM = mbsdManager->GetMBSDListM();
  fNSPState = mbsdManager->GetNSPState();
  fNParticle = mbsdManager->GetNParticle();
} // end of the constructor

TAHamiltonian::~TAHamiltonian(){} // end of the destructor

unsigned long TAHamiltonian::GetNBasis() const{
  return fMBSDListM->GetNBasis();
} // end of member function GetNBasis

/// calculate H(rr,cc) \param r: row, c: column
double TAHamiltonian::MatrixElement(int rr, int cc){
  return MatrixElement1N(rr, cc) + MatrixElement2N(rr, cc) +
    MatrixElement3N(rr, cc);
} // end of member function MatrixElement

void TAHamiltonian::DotProduct(TASparseVec &v, TASparseVec &r){
  const int n = GetNBasis();
  r.Clear(); // release the
  double t, e1, e2;
  for(int i = 0; i < n; i++){
    t = 0.;
    for(int j = 0; j < n; j++)
      if(!(e1 = MatrixElement(i,j)) && (!(e2 = v[j]))) t += e1 * e2;
    r.Fill(i,t);
  } // end for over i
} // end member function DotProduct
