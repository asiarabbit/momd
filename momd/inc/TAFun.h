/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAFun.h
  \class TAFun
  \brief To define a functor so as to be used in template functions that require
  a functor to be passed as an argument.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/25
  \date Last modified: 2020/07/25 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAFun_h_
#define _TAFun_h_

template<class T>
class TAFun{
public:
  TAFun(){}
  ~TAFun(){}

  virtual T operator()(double x) const = 0;
};

#endif
