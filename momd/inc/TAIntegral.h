/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAIntegral.h
  \class TAIntegral
  \brief This class implements general numerical integral. Input and results are
  all stored in arrays. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAIntegral_h_
#define _TAIntegral_h_

class TAIntegral{
public:
  TAIntegral(){}
  virtual ~TAIntegral(){}

  double Simpson(int n, const double *x, const double *f);
  // TAComplex Simpson(int n, const double *x, const TAComplex *f);
};

#endif
