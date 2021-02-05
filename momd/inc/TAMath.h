/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMath.h
  \class TAMath
  \brief Math class, to provide some general math methods.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/10/05 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMath_h_
#define _TAMath_h_

#include <cmath>
#include "TAMatrix.h"

class TAMath{
public:
  TAMath(){}
  virtual ~TAMath(){}

  // physical functions
  static constexpr double Pi(){ return 3.14159265358979323846; }
  /// golden cut ratio
	static constexpr double Alpha(){ return 0.61803398874989484820458683436565; }
  /// rad per degree
	static constexpr double DEGREE(){ return 0.01745329251994329547; }
	static constexpr double Sqrt3(){ return 1.73205080756887719318; }
  static constexpr double uMeV(){ return 931.494061; } ///< u (C12=12.u) in MeV
  static constexpr double hbarc(){ return 197.32705; } ///< \retval hbar*c in MeV*fm
  static constexpr double FineStructureConstant(){ return 7.29735307964482e-3; } ///< 1./137.0359895
  static constexpr double Pion0Mass(){ return 134.9770; } ///< \retval pion0 mass in MeV
  static constexpr double PionPMMass(){ return 139.5706; } ///< \retval pion+- mass in MeV
  static constexpr double c0(){ return 299792458.; } ///< \reteval light speed in m/s
  /// from Ek per u to velocity for a generic nucleon of mass uMeV
  static double EkPerUToBeta(double pEk){
    const double r = 1. + uMeV()/pEk;
    return sqrt(2.*r-1.)/r;
  }
  /// from Ek per u to momentum for a generic nucleon of mass uMeV
  static double EkPerUToMomentum(double pEk){
    const double b = EkPerUToBeta(pEk), g = 1./sqrt(1.-b*b);
    return b*g*uMeV();
  }
  static double WoodsSaxon(double r, double a0, double r0){ return 1./(1.+exp((r-r0)/a0)); }

  /// sign of a number
  static double sign(double c){
    if(c >= 0.) return 1.;
    return -1.;
  }
  /// sign of a number
  static double sign(int c){ return sign(double(c)); }
  /// \retval (-)^m
  static double minus(int n){ return n & 1 ? -1 : 1; }
  static double sqr(double a){ return a == 0. ? 0. : a*a; }
  static double innerProduct(int n, const double *p0, const double *p1){
    double prod = 0.;
    for(int i = n; i--;) prod += p0[i]*p1[i];
    return prod;
  }
  static double norm(int n, const double *p){ return sqrt(innerProduct(n, p, p)); }
  /// \retval computes (a^2+b^2)^0.5 without desctructive underflow or overflow
  static double norm(double a, double b){
    a = fabs(a); b = fabs(b);
    if(a > b) return a*sqrt(1.+sqr(b/a));
    else return 0. == b ? 0. : b*sqrt(1.+sqr(a/b));
  }
  static double infNorm(int n, const double *f){
    double m = 0.;
    for(int i = n; i--;) if(fabs(f[i]) > m) m = fabs(f[i]);
    return m;
  }
  // generic programming :)
  /// \retval: sum of square of each one in the parameter list
  static double Sum2(){ return 0.; }
  template <typename T, typename ...Args> // typename preferred, class deprecated
  static T Sum2(const T &v, const Args &...args){
    return v*v + Sum2(args...);
  }
  template <typename T>
  static bool Within(const T &v, const T &A, const T &B){
  if(v >= A && v <= B) return true;
  if(v >= B && v <= A) return true;
    return false;
  }

  /// \retval calculate combination number
  /// returns the value ln[Gamma(xx)] for xx > 0.
  static double gammaln(double xx);
  // returns ln(n!)
  static double factln(int n);
  /// \retval n!
  static int Factorial(int x);
  /// \param n, m: n should not be smaller than m
  static int Binomial(int n, int m);
  /// beta function, Binomial coefficients with float arguments
  static double Beta(double z, double w);
  /// \retval n!!
  static int BiFactorial(int n);
  /// the incomplete gamma function gammap, gammpa(a,x)=P(a,x)=gamma(a,x)/Gamma(a)
  static double gammap(double a, double x);
  /// the incomplete gamma function gammap, gammpa(a,x)=P(a,x)=gamma(a,x)/Gamma(a)
  static double gammaq(double a, double x);
  /// returns the incomplete gamma function P(a,x) evaluated by its series
  /// representation. Optionally returns ln[Gamma(a)] as gln.
  static double gser(double a, double x, double *gln = nullptr);
  /// returns the incomplete gamma function P(a,x) evaluated by its continued fraction
  /// representation. Optionally returns ln[Gamma(a)] as gln.
  static double gcf(double a, double x, double *gln = nullptr);

  /////// Linear algebra /////
  /// solve linear equation set a*x=b by Gauss-Jordan elimination
  /// \param a is a na*na matrix, and b is a na*nb matrix
  /// a is inversed inplace, and the solution is stored in b
  static void GaussJordan(matrix &a, int na, matrix &b, int nb);
  /// \retval calls GuassJordan(a,na,b,1)
  static void GaussJordan(matrix &a, int na, double *b);
  /// matrix inversion using Gauss-Jordan elimination: just inverse matrix a
  static void GaussJordan(matrix &a, int na){ matrix b(1, 1); GaussJordan(a, na, b, 0); }

  /// Given a matrix a[0...n-1][0...n-1], this routine replaces it by the LU decomposition
  /// of a rowwise permutation of itself. index[0...n-1] is an output vector that
  /// that records the row permutation effected by the partial pivoting; d is output
  /// as +-1 depending on whether the number of row interchanges was even or odd, respectively.
  /// This routine is used in combination with LUBackSubstitution to solve linear equations
  /// or invert a matrix
  static void LUDecomposition(matrix &a, int n, int *index, double &d);

  /// Solves the equation set ax=b. Here a[0...n-1][0...n-1] is input, not as the
  /// matrix A but rather as its LU decomposition, determined by the routine LUDecomposition.
  /// index[0...n-1] is input as the permutation vector returned by LUDecomposition.
  /// b[0..n-1] is input as the right-hand side and updated with the solution inplace.
  /// This routine takes into account the possibility that b will begin with many
  /// zeros, so it is efficient for use in matrix inversion
  static void LUBackSubstitution(const matrix &a, int n, const int *index, double *b);
  /// solve n-dim linar equation set ax=b using LU decomposition
  static void LUSolve(matrix &a, int na, double *b);
  static void LUSolve(matrix &a, int na, matrix &b, int nb);
  /// inverse n-dim matrix a using LU decomposition: AX=E
  static void LUInverse(matrix &a, int n, matrix &e);
  static void LUInverse(matrix &a, int n);
  /// calculate determinent of a square matrix. NOTE that this will change a
  static double Det(matrix &a, int n);

  /// expand in storage the covariance matrix covar, so as to take into account
  /// parameters that are being held fixed, which would be assigned zero covariances
  /// \param covar should be a ma*ma 2-D array
  /// \param nFit is the count of fitted parameters (not held fixed), -1: not assigned
  static void FillCovar(matrix &covar, int ma, bool *isFit, int nFit = -1);
};

#endif
