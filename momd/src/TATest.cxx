/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TATest.h
  \class TATest
  \brief Just to accommodate catch2 to facilitate unit test.
  It takes forever to compile this class. So normally leave this file alone.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/09/07
  \date Last modified: 2020/09/07 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include "TATest.h"

int TATest::test(int argc, char *argv[]){
  return Catch::Session().run(argc, argv);
} // end of member function test
