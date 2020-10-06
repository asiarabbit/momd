/**
	SUNNY Project, Anyang Normal University, IMP-CAS
	\file test.cxx
	\brief Just for general unit test
  \date Last modified: 2020/02/25
	\date Created: 2020/02/25
	\copyright 2020, SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include "TAMathFCI.h"

int main(){
  const int n = 2;
  TAMatrix2D ma(n, n); ma = {3, 2, 4, 5};
  TAMatrix2D v(n); v = {1., 1.};
  ma.Print();
  v.Print();

  TAMathFCI::EigenPower(ma, v);
}
