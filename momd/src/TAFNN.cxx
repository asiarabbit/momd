/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFNN.cxx
  \class TAFNN
  \brief This class represents nucleon-nucleon scattering amplitude, extracted
  from fitting of pp, pn and nn scattering data over a wide range of energies.
  fNN(q)=kNN/4pi*sigmaNN(i+alphaNN)*exp(-betaNN*q^2)
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TAFNN.h"
#include "TAComplex.h"
#include "TAInterplate.h"

/// \param kNN: wave number of n-n scattering
TAComplex TAFNN::GetFNN(double kNN){

}
