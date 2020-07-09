/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMOMDIS_M.cxx
  \class TAMOMDIS_M
  \brief Calculate core momentum distribution with angular momentum component m
  specified. This is a class to assist class TAMOMDIS.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/07/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TAMOMDIS_M.h"

TAMOMDIS_M::TAMOMDIS_M(TAMOMDIS *mom) : fMOM(mom){

}

TAMOMDIS_M::~TAMOMDIS_M(){}

/// returning c.s. for m, mom is assigned to momStr
double TAMOMDIS_M::ParallelStr(int m, double *momStr){

}
/// the same as ParallelStr, but for diffraction dissociation
double TAMOMDIS_M::ParallelDiff(int m, double *momDiff){

}
