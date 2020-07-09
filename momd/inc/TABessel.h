/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TABessel.h
  \class TABessel
  \brief Numerical result of Bessel functions, stored in arrarys. Note that this
  class is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TABessel_h_
#define _TABessel_h_

class TABessel{
public:
  TABessel(){}
  virtual ~TABessel(){}

  // Bessel function of integer order
  static double BesselJ(int n, double x); // integer order
  static double BesselJ(double a, double x); // half integer order
  // modified Bessel function: I_a(x) = i^{-a}*J_a(ix)
  static double BesselI(int n, double x); // integer order
  static double BesselI(double a, double x); // half integer order
};

#endif
