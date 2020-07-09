/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAComplex.h
  \class TAComplex
  \brief A datatype for a complex double number.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAComplex_h_
#define _TAComplex_h_

#include <complex>

class TAComplex : public complex<double>{
public:
  TAComplex(double real, double imag);
  virtual ~TAComplex(){}
};

#endif
