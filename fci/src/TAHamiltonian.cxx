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
#include "TASingleParticleState.h"
#include "TASingleParticleStateManager.h"

static const double DUMB = -9999.;

TAHamiltonian::TAHamiltonian(const string &configFile)
    : fMBSDListM(0), fMatrix(0), fNSPState(0), fNMBSD(0){
  // prepare the single-particle basis set //
  // prepare the basis of the representation //
  TAManyBodySDManager *mbsdManager = TAManyBodySDManager::Instance();
  mbsdManager->LoadConfigFile(configFile);
  mbsdManager->MSchemeGo(); // generate many-body basis
  fMBSDListM = mbsdManager->GetMBSDListM();
  fNSPState = mbsdManager->GetNSPState();
  fNParticle = mbsdManager->GetNParticle();
  fNMBSD = fMBSDListM->GetNBasis();
} // end of the constructor

TAHamiltonian::~TAHamiltonian(){
  if(fMatrix){ delete fMatrix; fMatrix = nullptr; }
} // end of the destructor

/// \retval calculate and return the matrix form of the hamiltonian
matrix* TAHamiltonian::Matrix(){
  if(fMatrix && !fMatrix->IsEmpty()){
    return fMatrix;
  }
  if(!fMBSDListM)
    TAException::Error("TAHamiltonian",
      "Matrix: The many-body basis not assigned.");
  if(!fNSPState)
    TAException::Error("TAHamiltonian", "Matrix: fNSPState not assigned.");
  if(!fNMBSD)
    TAException::Error("TAHamiltonian", "Matrix: fNMBSD is 0. Not assigned?");

  // loop to generate each matrix element for the hamiltonian //
  if(fMatrix){ delete fMatrix; fMatrix = nullptr; }
  fMatrix = new matrix(fNMBSD, fNMBSD); // allot memery to a nxn matrix
  // initialize to a specific initial value //
  for(int i = fNMBSD; i--;) for(int j = fNMBSD; j--;) (*fMatrix)[i][j] = DUMB;

  // H is symmetric - alculate the lower triangle, and then copy back //
  for(int rr = 0; rr < fNMBSD; rr++) for(int cc = 0; cc <= rr; cc++){
    MatrixElement(rr, cc); // assign matrix element H[i][j]
    (*fMatrix)[cc][rr] = (*fMatrix)[rr][cc];
  } // end for

  return fMatrix;
} // end of the member function Matrix

/// assign the matrix element (*fMatrix)[i][j]
/// \param r: row, c: column
void TAHamiltonian::MatrixElement(int rr, int cc){
  // so that the matrix is symmetric //
  if((*fMatrix)[rr][cc] != DUMB)
    TAException::Warn("TAHamiltonian", "H element (%d, %d) assgiend already assigned.", rr, cc);
  (*fMatrix)[rr][cc] = MatrixElement1N(rr, cc) + MatrixElement2N(rr, cc) +
     MatrixElement3N(rr, cc);
} // end of member function MatrixElement
