/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TALegendre.cxx
  \class TALegendre
  \brief To calculate associate Legendre polynomials for given l and m. Results
  are stored in arrays. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/09/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cmath>
#include <catch2/catch.hpp>
#include "TALegendre.h"
#include "TAMath.h"
#include "TAException.h"


double TALegendre::Legendre(int l, int m, double x){
  if(m >= 0) return LegendreM(l, m, x);
  return LegendreM(l, -m, x) *
    exp(TAMath::factln(l+m)-TAMath::factln(l-m))*TAMath::minus(m);
} // end of member function Legendre
double TALegendre::LegendreM(int l, int m, double x){
  if(l < 0.) TAException::Error("TALegendre",
    "Legendre(...): input l is minus.");
  if(fabs(x) > 1.) TAException::Error("TALegendre",
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
} // end of member function LegendreM

TEST_CASE("Legendre Polynomials", "[legendre]"){
  CHECK(TALegendre::Legendre(0,0,1.) == 1.);
  CHECK(TALegendre::Legendre(1,0,1.) == 1.);
  CHECK(TALegendre::Legendre(7,0,0.5) == Approx(0.22314453125).epsilon(1.e-10));
  CHECK(TALegendre::Legendre(5,2,0.5) == Approx(-4.921875).epsilon(1.e-10));
  CHECK(TALegendre::Legendre(5,-2,0.5) == Approx(-0.005859375).epsilon(1.e-10));
  CHECK(TALegendre::Legendre(5,-2,-0.5) == Approx(0.005859375).epsilon(1.e-10));
  CHECK(TALegendre::Legendre(5,-4,0.5) == Approx(0.000732421875).epsilon(1.e-10));
  CHECK(TALegendre::Legendre(5,-3,0.3) == Approx(-0.0004295210623).epsilon(1.e-10));
}
