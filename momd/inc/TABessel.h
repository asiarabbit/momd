/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TABessel.h
  \class TABessel
  \brief Numerical results of Bessel and modified Bessel functions of integer order.
  Note that this class is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/09/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TABessel_h_
#define _TABessel_h_

class TABessel{
public:
  TABessel(){}
  virtual ~TABessel(){}

  /// the implementation here are directly taken from: Numerical Receipes in C
  /// Bessel function of integer order
  /// BesselJ0,1 are computed first. The other orders are calculated using the
  /// the recurrence relations.
  static double BesselJ0(double x);
  static double BesselY0(double x);
  static double BesselJ1(double x);
  static double BesselY1(double x);
  static double BesselJ(int n, double x); // integer order
  static double BesselY(int n, double x); // integer order

  // modified Bessel functions
  static double BesselI0(double x);
  static double BesselK0(double x);
  static double BesselI1(double x);
  static double BesselK1(double x);
  static double BesselI(int n, double x);
  static double BesselK(int n, double x);
};

#endif
