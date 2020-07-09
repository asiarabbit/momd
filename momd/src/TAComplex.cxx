/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAComplex.cxx
  \class TAComplex
  \brief A datatype for a complex double number.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TAComplex.h"

using std::complex;

class TAComplex : public complex<double>{
public:
  TAComplex(double real, double imag) : complex<double>(real, imag){}
  virtual ~TAComplex(){}
};

#endif
