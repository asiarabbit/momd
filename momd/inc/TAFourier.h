/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFourier.h
  \class TAFourier
  \brief This class implements Fourier transfrom of f(r) to f(q). Input and
  results are all stored in arrays. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAFourier_h_
#define _TAFourier_h_

class TAComplex;

class TAFourier{
public:
  TAFourier(){}
  virtual ~TAFourier(){}

  // 3D Cartesian space Fourier transform in transverse (xy) plane
  // input func is fr; output is fq
  static void Fourier2D(const double *fr, TAComplex *fq);
};

#endif
