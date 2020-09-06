/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TALegendre.cxx
  \class TALegendre
  \brief To calculate associate Legendre polynomials for given l and m. Results
  are stored in arrays. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/09/03 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TALegendre.h"
#include "TAMath.h"
#include "TAException.h"

double TALegendre::Legendre(int l, int m, double x){
  if(l < 0.) TAException::Error("TALegend",
    "Legendre(...): input l is minus.");
  if(fabs(x) > 1.) TAException::Error("TALegend",
    "Legendre(...): input x is not within [-1,1].");
  if(l < m) return 0.;

  // Pmm =(-)^m*(2m-1)!!(1-x^2)^(m/2)
  double pmm = TAMath::minus(m) * TAMath::BiFactorial(2.*m-1.)
    * pow((1.+x)*(1.-x), m/2.);
  if(m == l) return pmm;
  double pmm1 = x*(2.*m+1.)*pmm; // P_(m+1)^m
  if(m+1 == l) return pmm1;

  // recurrence (i+2-m)P(i+2)m=x(2li+3)P(i+1)m-(i+m+1)Pim
  double pmm2 = 0.;
  for(int i = m + 2; i <= l; i++){
    pmm2 = (x*(2*i-1.)*pmm1-(i+m-1.)*pmm)/(i-m);
    pmm = pmm1; pmm1 = pmm2;
  } // end for over i
  return pmm2;
} // end of member function Legendre(...)
