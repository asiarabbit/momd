/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMath.h
  \class TAMath
  \brief Math class, to provide some general math methods. Note that this is a
  tool class, so it is defined to be a static class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/07/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TAMath.h"

/// sign of a number
double TAMath::sign(double c){
  if(c >= 0.) return 1.
  return -1.;
} // end of member function sign
