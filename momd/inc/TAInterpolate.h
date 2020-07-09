/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAInterpolate.h
  \class TAInterpolate
  \brief This is a general interpolation class. The principle method used is
  Lagrange interpolation method. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/07/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAInterpolate_h_
#define _TAInterpolate_h_

class TAInterplate{
public:
  TAInterpolate(){}
  virtual ~TAInterpolate(){}

  /// \param len: length of array x and y; \param order: order of the polynomial
  static double Lagrange(double *x, double *y, int len, int order = 5);
};
