/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TATest.h
  \class TATest
  \brief Just to accommodate catch2 to facilitate unit test
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/09/07
  \date Last modified: 2020/09/07 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TATest_h_
#define _TATest_h_

class TATest{
public:
  TATest(){}
  virtual ~TATest(){}

  static int test(int argc, char *argv[]);
};

#endif
