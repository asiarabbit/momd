/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TAFCI.h
  \class TAFCI
  \brief Global class for full configuration interaction calculations, including
  but not limited to read-in of user-input single-particle states, generation of
  many-body basis, establishment and diagonalization of the Hamiltonian matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/02/12
  \date Last modified: 2020/02/12 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAFCI_h_
#define _TAFCI_h_

class TAHamiltonian;

class TAFCI{
public:
  virtual ~TAFCI();
  static TAFCI *Instance();

  void Go(); // initiate the FCI calculation

private:
  TAFCI();

  static TAFCI *kInstance;
  TAHamiltonian *fHamiltonian;
};

#endif
