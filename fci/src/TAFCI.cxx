/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TAFCI.cxx
  \class TAFCI
  \brief Global class for full configuration interaction calculations, including
  but not limited to read-in of user-input single-particle states, generation of
  many-body basis, establishment and diagonalization of the Hamiltonian matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/02/12
  \date Last modified: 2020/10/07 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include "TAFCI.h"
#include "TAHamiltonian.h"
#include "TAMath.h"

TAFCI *TAFCI::kInstance = nullptr;

TAFCI::TAFCI() : fHamiltonian(nullptr){}

TAFCI::~TAFCI(){}

TAFCI *TAFCI::Instance(){
  if(!kInstance) kInstance = new TAFCI();
  return kInstance;
} // end of member function Instance

void TAFCI::Go(){
//  if(!fHamiltonian) fHamiltonian = TAHamiltonian::Instance();

//  fHamiltonian->InitializeCoefficient(); // DEBUG
//  matrix &H = fHamiltonian->Matrix();
//  H.Print();
//  H.PrintInC();
//  getchar(); // DEBUG

  const int n = 1; // H.nc();
  matrix H(n, n);
  matrix v(n);
  for(int i = n; i--;) v[i][0] = 1.;

  // solve the dominant eigenvalue and eigenvector using power method //

  // solve all of the eigenvalues using QR method //
  // firstly, let's implement QR factorization for H
//  matrix Q(n, n), R(n, n);
//  TAMath::QR(H, Q, R);

  // QR method to calculate all the eigenvalues of a matrix
  matrix vv(n);
//  return;
  // Jacobi method to calculate all the eigenvalues and eigenvectors of a matrix
  matrix P(n, n);
//  for(int i = H.nr(); i--;) H[i][i] -= 50.;
//  TAMath::EigenPower(H, v);
//  TAMath::EigenQR(H, vv);
//  TAMath::EigenJacobi(H, P, vv);
return;
  // TAMath::EigenDavidson(3, H, P, vv);

} // end of member function Go
