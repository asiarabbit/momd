/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TABessel.cxx
  \class TABessel
  \brief Numerical result of Bessel functions, stored in arrarys. Note that this
  class is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/07/19 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TABessel.h"
#include "TAMath.h"
#include "TAException.h"

// Bessel function: 0-order the first kind
// Numerical Receipes in C: p232 (ISBN: 0-521-43108-5)
// x: [0,8]: rational functions of x
// x: [8: +\infty]: Jn(x)=sqrt(2/(PI*x))[Pn(8/x)cos(Xn)-Qn(8/x)sin(Xn)]
// Xn: x - (2n+1)PI/4
double TABessel::BesselJ0(double x){
  x = fabs(x);

  if(x < 8.){
    double y = x*x;
    // direct rational function fit
    double y1 = 57568490574.+y*(-13362590354.+y*(651619640.7
      +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    double y2 = -0.1562499995e-1+y*(0.1430488765e-3
      +y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)));
    return y1 / y2;
  } // end if

  // x >= 8.
  double z = 8./x, y = z*z;
  double xn = x - 0.785398163397448; // x - PI/4
  double pn = 1.+y*(-0.1098628627e-2+y*(0.2734510407e-4
    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
  double qn = -0.1562499995e-1+y*(0.1430488765e-3
    +y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)));
  qn *= z;
  // sqrt(2/Pi/x)
  return sqrt(0.6366197723676/x) * (pn*cos(xn) - qn*sin(xn));
} // end of member function BesselJ0

// order-0 the second kind: 0-th Neumann function
double TABessel::BesselY0(double x){
  if(x <= 0.)
    TAException::Error("TABessel",
      "BesselY0(): Only defined for positive x.\n");

  if(x < 8.0){
    double y = x*x;
    double y1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
      +y*(10879881.29+y*(-86327.92757+y*228.4622733))));
    double y2 = 40076544269.0+y*(745249964.8+y*(7189466.438
      +y*(47447.26470+y*(226.1030244+y*1.0))));
    return y1/y2 + 0.6366197723676*BesselJ0(x)*log(x);
  } // end if

  // x >= 8.
  double z = 8./x, y = z*z;
  double xn = x - 0.785398163397448; // x - PI/4
  double pn = 1.+y*(-0.1098628627e-2+y*(0.2734510407e-4
    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
  double qn = -0.1562499995e-1+y*(0.1430488765e-3
    +y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)));
  qn *= z;
  // sqrt(2/Pi/x)
  return sqrt(0.6366197723676/x) * (pn*sin(xn) + qn*cos(xn));
} // end of member function BesselY0

double TABessel::BesselJ1(double x){
  double sign = TAMath::sign(x);
  x = sign;

  if(x < 8.){
    double y = x*x;
    // direct rational function fit
    double y1 =x*(72362614232.+y*(-7895059235.+y*(242396853.1
      +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    double y2 = 144725228442.+y*(2300535178.+y*(18583304.74
      +y*(99447.43394+y*(376.9991397+y*1.))));
    return sign * y1 / y2;
  } // end if

  // x >= 8.
  double z = 8./x, y = z*z;
  double xn = x - 2.356194490192345; // x - (2+1)PI/4
  double pn = 1.+y*(0.183105e-2+y*(-0.3516396496e-4
    +y*(0.2457520174e-5+y*(-0.240337019e-6))));
  double qn = 0.04687499995+y*(-0.2002690873e-3
    +y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
  qn *= z;
  // sqrt(2/Pi/x)
  return sign * sqrt(0.6366197723676/x) * (pn*cos(xn) - qn*sin(xn));
} // end of member function BesselJ1

// order-0 the second kind: 0-th Neumann function
double TABessel::BesselY1(double x){
  if(x <= 0.)
    TAException::Error("TABessel",
      "BesselY1(): Only defined for positive x.\n");

  if(x < 8.0){
    double y = x*x;
    double y1 = x*(-0.4900604943e13+y*(0.1275274390e13
      +y*(-0.5153438139e11+y*(0.7349264551e9
      +y*(-0.4237922726e7+y*0.8511937935e4)))));
    double y2 = 0.2499580570e14+y*(0.4244419664e12
      +y*(0.3733650367e10+y*(0.2245904002e8
      +y*(0.1020426050e6+y*(0.3549632885e3+y)))));
    return y1/y2 + 0.6366197723676*(BesselJ1(x)*log(x) - 1./x);
  } // end if

  // x >= 8.
  double z = 8./x, y = z*z;
  double xn = x - 2.356194490192345; // x - (2+1)PI/4
  double pn = 1.+y*(0.183105e-2+y*(-0.3516396496e-4
    +y*(0.2457520174e-5+y*(-0.240337019e-6))));
  double qn = 0.04687499995+y*(-0.2002690873e-3
    +y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
  qn *= z;
  // sqrt(2/Pi/x)
  return sqrt(0.6366197723676/x) * (pn*sin(xn) + qn*cos(xn));
} // end of member function BesselY1

double TABessel::BesselY(int n, double x){
  if(x <= 0.)
    TAException::Error("TABessel",
      "BesselY(): Only defined for positive x.\n");
  if(0 == n) return BesselY0(x);
  if(1 == n) return BesselY1(x);

  // calculate Yn using Y(n+1)=2n/x*Yn-Y(n-1)
  double byp, bym = BesselY0(x), by = BesselY1(x);
  double tox = 2./x; // two over x
  for(int i = 1; i < n; i++){
    byp = i*tox * by - bym;
    bym = by; by = byp;
  } // end for over i
  return by;
} // end of the member function BesselY

// calculate BesselJ using upward recurrence for x > n
// and downward recurrence for x <= n
// ACC: the number of significant figures of accuracy.
static const double ACC = 40.;
// BIGNO, BIGNI: renormalization factor to prevent overflows
static const double BIGNO = 1.e10, BIGNI = 1.e-10;
double TABessel::BesselJ(int n, double x){
  if(0 == n) return BesselJ0(x);
  if(1 == n) return BesselJ1(x);

  double sign = TAMath::sign(x);
  x = fabs(x);
  if(0. == x) return 0.;

  double tox = 2./x; // two over x
  double bjp, bjm = BesselJ0(x), bj = BesselJ1(x);
  if(x > double(n)){
    for(int i = 1; i < n; i++){
      bjp = i*tox * bj - bjm;
      bjm = bj; bj = bjp;
    } // end for over i
    return bj;
  } // end if

  // downward recurrence from an even nmax
  int nmax = 2.*( (n+int(sqrt(ACC*n)))/2. );
  double sum = 0., ans = 0.;
  // 1 = J0 + 2J2 + 2J4 + 2J6 + ...
  bool isAdd = false; // if put J into addition
  bjp = 0.; bj = 1.;
  for(int i = nmax; i > 0.; i--){
    bjm = i*tox *bj - bjp;
    bjp = bj; bj = bjm;
    if(fabs(bj) > BIGNO){ // renormalize to prevent overflows
      bj *= BIGNI; bjp *= BIGNI;
      sum *= BIGNI; ans *= BIGNI;
    } // end if
    if(isAdd) sum += bj;
    isAdd = !isAdd;
    if(i == n) ans = bj;
  } // end for over i
  sum = 2. * sum - bj; // J0 + 2J2 + 2J4 + 2J6 + ...
  ans /= sum; // normalize Jn
  return n & 1 ? sign * ans : ans;
} // end of the member function BesselJ

// modified Bessel function: I_n(x) = i^{-n}*J_n(ix)
// integer order
double TABessel::BesselI0(double x){
  x = fabs(x);

  double y = x / 3.75;
  if(x < 3.75){ // polynomial fit
    y *= y;
    return 1.+y*(3.5156229+y*(3.0899424+y*(1.2067492
      +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  }
  y = 3.75 / x;
  return
    (exp(x)/sqrt(x))*(0.39894228+y*(0.1328592e-1
    +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
    +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
    +y*0.392377e-2))))))));
} // end of member function BesselI0

// modified Bessel function: I_n(x) = i^{-n}*J_n(ix)
// integer order
double TABessel::BesselK0(double x){
  if(x <= 0.)
    TAException::Error("TABessel",
      "BesselK0(): Only defined for positive x.\n");

  double y = 0.;
  if(x <= 2.){ // polynomial fit
    y = x*x / 4.;
    return (-log(x/2.)*BesselI0(x))+(-0.57721566+y*(0.42278420
      +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
      +y*(0.10750e-3+y*0.74e-5))))));
  }
  y = 2. / x;
  return
    (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
    +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
    +y*(-0.251540e-2+y*0.53208e-3))))));
} // end of member function BesselK0

double TABessel::BesselI1(double x){
  double sign = TAMath::sign(x);
  x = fabs(x);

  double y = x / 3.75;
  if(x < 3.75){ // polynomial fit
    y *= y;
    return x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
      +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  }
  y = 3.75 / x;
  double y1 = 0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
  y1 = 0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
    +y*(0.163801e-2+y*(-0.1031555e-1+y*y1))));
  return y1 * (exp(x)/sqrt(x)) * sign;
} // end of member function BesselI1

// modified Bessel function: I_n(x) = i^{-n}*J_n(ix)
// integer order
double TABessel::BesselK1(double x){
  if(x <= 0.)
    TAException::Error("TABessel",
      "BesselK1(): Only defined for positive x.\n");

  double y = 0.;
  if(x <= 2.){ // polynomial fit
    y = x*x / 4.;
    return log(x/2.)*BesselI1(x)+(1./x)*(1.+y*(0.15443144
      +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
      +y*(-0.110404e-2-y*0.4686e-4))))));
  }
  y = 2. / x;
  return
    (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
    +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
    +y*(0.325614e-2-y*0.68245e-3))))));
} // end of member function BesselK1

// modified Bessel function of the 2nd kind
// upward recurrence for all x
double TABessel::BesselK(int n, double x){
  if(x <= 0.)
    TAException::Error("TABessel",

      "BesselK(): Only defined for positive x.\n");
  if(0 == n) return BesselK0(x);
  if(1 == n) return BesselK1(x);
  // calculate Kn using K(n+1)=2n/x*Kn-K(n-1)
  double bkp, bkm = BesselK0(x), bk = BesselK1(x);
  double tox = 2./x; // two over x
  for(int i = 1; i < n; i++){
    bkp = i*tox * bk + bkm;
    bkm = bk; bk = bkp;
  } // end for over i
  return bk;
} // end of the member function BesselK

// modified Bessel function of the first kind
// upward recurrence for x > n
// and downward recurrence for x <= n
// ACC: the number of significant figures of accuracy.
double TABessel::BesselI(int n, double x){
  static const double BIGNO = 1.e10, BIGNI = 1.e-10;

  if(0 == n) return BesselI0(x);
  if(1 == n) return BesselI1(x);

  double sign = TAMath::sign(x);
  x = fabs(x);
  if(0. == x) return 0.;

  double tox = 2./x; // two over x
  double bip = 0., bi = 1., bim, ans;
  // downward recurrence from an even nmax
  int nmax = 2.*( n+int(sqrt(ACC*n)) );
  // 1 = I0 - 2J2 + 2J4 - 2J6 + ...
  for(int i = nmax; i > 0.; i--){
    bim = bip - i*tox *bi;
    bip = bi; bi = bim;
    if(fabs(bi) > BIGNO){ // renormalize to prevent overflows
      bi *= BIGNI;
      bip *= BIGNI;
      ans *= BIGNI;
    } // end if
    if(i == n) ans = bi;
  } // end for over i
  ans *= BesselI0(0) / bi; // normalize Jn
  return n & 1 ? sign * ans : ans;
} // end of the member function BesselI
