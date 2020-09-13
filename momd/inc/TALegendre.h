/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TALegendre.h
  \class TALegendre
  \brief To calculate associate Legendre polynomials for given l and m. Results
  are stored in arrays. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/19 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TALegendre_h_
#define _TALegendre_h_

class TALegendre{
public:
  TALegendre(){}
  virtual ~TALegendre(){}

  static double Legendre(int l, int m, double x); ///< for arbitrary integer m
  static double LegendreM(int l, int m, double x); ///< for m >= 0
};

#endif
