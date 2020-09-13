/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMath.h
  \class TAMath
  \brief Math class, to provide some general math methods. Note that this is a
  tool class, so it is defined to be a static class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/09/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cmath>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <utility>
#include <catch2/catch.hpp>
#include "TAMath.h"
#include "TAException.h"

using std::swap;
using std::max;

/// sign of a number
double TAMath::sign(double c){
  if(c >= 0.) return 1.;
  return -1.;
} // end of member function sign

/// returns the value ln[Gamma(xx)] for xx > 0.
double TAMath::gammaln(double xx){
  if(xx <= 0.) TAException::Error("TAMath", "gammaln: x > 0. is mandatory.");
  // internal arithmetic will be done in double precision, a nicety that you can
  // omit if five-figure accuracy is good enough
  static const double c[6] = {76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
  double tmp = xx + 5.5, x = xx; // (z+gamma+1/2), gamma = 5, N = 6
  tmp -= (x+0.5)*log(tmp);
  double ser = 1.000000000190015; // the additive series
  for(int i = 0; i < 6; i++) ser += c[i]/++x;
  return -tmp+log(2.5066282746310005*ser/xx); // sqrt(2pi)
} // end of member function gammaln

// returns ln(n!)
double TAMath::factln(int n){
  static double a[101]; // a static array is automatically initialized to zero

  if(n < 0) TAException::Error("TAMath", "factln: n is negative.");
  if(n <= 1) return 0.;
  if(n <= 100) return a[n] ? a[n] : (a[n]=gammaln(n+1.));
  return gammaln(n+1.); // out of range of table
} // end of member function factln

/// \retval n!
int TAMath::Factorial(int n){
  static int ntop = 4;
  static int a[33] = {1, 1, 2, 6, 24}; // fill in table only as required

  if(n < 0) TAException::Error("TAMath", "Factorial: n: %d is minus", n);
  // larger value than size of the table is required. Actually, this big a value
  // is going to overflow on many computers, but no harm in trying
  if(n > 32) return exp(gammaln(n+1.));
  while(ntop < n){
    ntop++;
    a[ntop] = a[ntop-1]*ntop;
  }
  return floor(0.5+a[n]); // clear off the roundoff error
} // end of member function Factorial

// returns Binomial coefficients
int TAMath::Binomial(int n, int m){
  if(m > n)
    TAException::Error("TAMath",
      "Binomial: n: %d is larger than m: %d!", n, m);

  // the floor function cleans up roundoff error for smaller values of n and m
  return floor(0.5+exp(factln(n)-factln(m)-factln(n-m)));
} // end of member function Binomial

double TAMath::Beta(double z, double w){
  if(z < 0. || w < 0.) TAException::Error("TAMath", "Beta: z or z is minus");
  return exp(gammaln(z)+gammaln(w)-gammaln(z+w));
} // end of member function Beta

/// \retval n!!
int TAMath::BiFactorial(int n){
  return n <= 1 ? 1 : n * BiFactorial(n-2);
} // end of member function BiFactorial

/// the incomplete gamma function gammap, gammap(a,x)=P(a,x)=gamma(a,x)/Gamma(a)
double TAMath::gammap(double a, double x){
  if(x < 0. || a <= 0.) TAException::Error("TAMath", "gammap: invalid arguements");

  if(x < (a+1.)) return gser(a, x); // series representation
  return 1. - gcf(a, x);
} // end of member function gammap
/// the incomplete gamma function gammaq, gammaq(a,x)=Q(a,x)=1-P(a,x)
double TAMath::gammaq(double a, double x){
  if(x < 0. || a <= 0.) TAException::Error("TAMath", "gammaq: invalid arguements");

  if(x < (a+1.)) return 1. - gser(a, x); // continued fraction representation
  return gcf(a, x);
} // end of member function gammaq

TEST_CASE("Incomplete Gamma function", "[gammap]"){
  // CHECK(TAMath::gammap(1., 0.5) == Approx(0.606531).epsilon(1.e-8));
  // CHECK(TAMath::gammaq(1., 0.5) == Approx(1.-0.606531).epsilon(1.e-8));
  CHECK(TAMath::gammap(1., 0.5) == Approx(0.3934693402873666).epsilon(1.e-8));
  CHECK(TAMath::gammap(1., 0.1) == Approx(0.09516258196404048).epsilon(1.e-8));
  CHECK(TAMath::gammaq(1., 3.5) == Approx(0.0301973834223185).epsilon(1.e-8));
  CHECK(TAMath::gammaq(1., 4.1) == Approx(0.016572675401761255).epsilon(1.e-8));
  CHECK(exp(TAMath::gammaln(0.5)) == Approx(sqrt(TAMath::Pi())).epsilon(1.e-8));
}

// returns the incomplete gamma function P(a,x) evaluated by its series
// representation. Optionally returns ln[Gamma(a)] as gln.
double TAMath::gser(double a, double x, double *gln){
  static const double ITMAX = 100; // maximum iteration times
  static const double EPS = 3.e-7; // fractional accuracy

  if(x < 0.) TAException::Error("TAMath", "gser: x < 0., x: %f", x);
  double ap = a, sum, item = sum = 1./a;
  for(int i = 0; i < ITMAX; i++){
    sum += (item *= x/++ap);
    if(fabs(item) < fabs(sum)*EPS)
      return sum * exp(-x + a*log(x) - (gln ? *gln=gammaln(a) : gammaln(a)));
  } // end for over i
  TAException::Error("TAMath", "gser: a too large, and IMAX too small.");
  return 0.; // never gets here
} // end of member function gamser

// returns the incomplete gamma function Q(a,x) evaluated by its continued
// fraction representation. Optionally returns ln[Gamma(a)] as gln.
double TAMath::gcf(double a, double x, double *gln){
  static const double ITMAX = 100; // maximum iteration times
  static const double EPS = 3.e-7; // factional accuracy
  // number near the smallest representable floating-point number
  static const double FPMIN = 1.e-30; // modified Lentz's algorithm

  // set up for evaluating continued fraction by modified Lentz's method with b0=0.
  // b1 = x+1.-a, c1 = \infty (since c0 = 0), d1 = 1/b1, fr1=b0+a1/b1=a1/b1=1/b1=d1
  double b = x+1.-a, c = 1./FPMIN, d = 1./b, fr = d;
  for(int i = 1; i < ITMAX; i++){
    double aip = -i*(i-a); // a_{i+1}
    b += 2.; d = aip*d+b; // b2, d2 = 1/(a2*d1+b2)
    if(fabs(d) < FPMIN) d = FPMIN;
    c = b+aip/c; // c2 = b2 + a2/c1
    if(fabs(c) < FPMIN) c = FPMIN;
    d = 1./d;
    double delta = c*d;
    fr *= delta; // fr_2 = fr_1*c2*d2
    if(fabs(delta-1.) < EPS)
      return exp(-x+a*log(x)-(gln ? *gln = gammaln(a) : gammaln(a)))*fr; // put factors in front
  } // end for over i
  TAException::Error("TAMath", "gcf: a too large, and IMAX too small.");
  return 0.; // never gets here
} // end of member function gcf

/////// Linear algebra /////
/// solve linear equation set by Gauss-Jordan elimination
/// \param a should be a ma*ma 2-D array
/// solve linear equation set a*x=b by Gauss-Jordan elimination
/// \param a is a na*na matrix, and b is a na*nb matrix
/// a is inversed inplace, and the solution is stored in b
/// full-pivoting is used (row & column)
void TAMath::GaussJordan(matrix &a, int na, matrix &b, int nb){
  int pivotRow[na], pivotCol[na], isPivoted[na]{}; // (row, col) of each pivot in stage i
  memset(pivotRow, -1, sizeof(pivotRow));
  memset(pivotCol, -1, sizeof(pivotCol));

  //// the main loop over the columns to be reduced
  for(int i = 0; i < na; i++){ // the i-stage of the algorithm
    double pivot = 0.;
    int irow, icol; // position of the pivot
    // search pivot globally (the one with the largest magnitude) //
    // pivot cannot be selected from the pivoted cols and rows
    // because it would mess up the built-up part of the unity matrix in matrix a
    for(int j = 0; j < na; j++) if(!isPivoted[j])
      for(int k = 0; k < na; k++) if(!isPivoted[k])
        if(fabs(a[j][k]) > pivot){ pivot = a[j][k]; irow = j; icol = k; }
    // the icol column has been pivoted
    if(fabs(pivot) <= 1E-100) TAException::Error("TAMath",
      "GaussJordan: Sigular coefficient matrix, pivot[%d]: %f", i, pivot);
    isPivoted[icol] = true;
    pivotRow[i] = irow; pivotCol[i] = icol;
    // interchange row irow and icol to put the pivot on the diagonal
    if(irow != icol){
      for(int j = 0; j < na; j++) swap(a[irow][j], a[icol][j]);
      for(int j = 0; j < nb; j++) swap(b[irow][j], b[icol][j]);
    } // end if
    // now we have the valid pivot, here the elimination begins //
    double pivotInv = 1./pivot;
    // rescale the icol row by 1./pivot
    a[icol][icol] = 1.; // to simulate the behavior of the identity matrix
    for(int j = 0; j < na; j++) a[icol][j] *= pivotInv;
    for(int j = 0; j < nb; j++) b[icol][j] *= pivotInv;
    for(int j = 0; j < na; j++){ // loop over rows
      if(j != icol){ // then subtract a[icol]*a[j][icol] from row-j of the augmented matrix
        double dum = a[j][icol];
        a[j][icol] = 0.; // so that the icol-column would end up, like a^-1, to be -a[j][icol]/pivot
        for(int k = 0; k < na; k++) a[j][k] -= dum*a[icol][k];
        for(int k = 0; k < nb; k++) b[j][k] -= dum*b[icol][k];
      } // end else
    } // end for over j
  } // end the i-stage
  // row interchange is not reflected on the resulted matrix, it is recovered by column interchange
  // (P2*P1*A)^-1 = B = A^-1*P1^-1*P2^-1 => A^-1=B*P2*P1
  for(int i = na; i--;) if(pivotRow[i] != pivotCol[i])
    for(int j = 0; j < na; j++) swap(a[j][pivotRow[i]], a[j][pivotCol[i]]);
} // end of member function GaussJordan-0
void TAMath::GaussJordan(matrix &a, int na, double *b){
  matrix bb(na, 1);
  for(int i = na; i--;) bb[i][0] = b[i];
  GaussJordan(a, na, bb, 1);
  for(int i = na; i--;) b[i] = bb[i][0];
} // end of member function GaussJordan-1

/// Given a matrix a[0...n-1][0...n-1], this routine replaces it by the LU decomposition
/// of a rowwise permutation of itself. index[0...n-1] is an output vector that
/// that records the row permutation effected by the partial pivoting; d is output
/// as +-1 depending on whether the number of row interchanges was even or odd, respectively.
/// This routine is used in combination with LUBackSubstitution to solve linear equations
/// or invert a matrix
void TAMath::LUDecomposition(matrix &a, int n, int *index, double &d){
  static const double TINY = 1.e-20;

  // loop over rows to get the implicit scaling information //
  d = 1.; // no row interchanges yet
  double big = 0., sum = 0., dum; // temporary variables
  double vv[n]; // vv: stores the implicit scaling of each row, 1/||row_i||_\infty
  for(int i = 0; i < n; i++){
    big = 0.;
    for(int j = 0; j < n; j++) big = max(big, fabs(a[i][j]));
    if(0. == big) TAException::Error("TAMath", "LUDecomposition: a is a Singular matrix.");
    vv[i] = 1./big; // save the scaling
  } // end for over i
  // loop over columns to implement Crout's method (LU decomposition) //
  for(int j = 0; j < n; j++){
    for(int i = 0; i < j; i++){ // solve U-betaij = aij-sum_(k=0,i-1){alphaik*betakj}
      sum = a[i][j];
      for(int k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    } // end for over i
    big = 0.; // initialize for the search for largest pivot element
    int imax; // imax: maximum beta in column imax
    for(int i = j; i < n; i++){
      sum = a[i][j];
      for(int k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if((dum = vv[i]*fabs(sum)) > big){ // search for the largest betajj
        big = dum; imax = i;
      } // end if
    } // end for over i
    if(j != imax){ // do we need to interchange rows?
      for(int k = 0; k < n; k++) swap(a[j][k], a[imax][k]);
      d = -d; vv[imax] = vv[j]; // and together with the parity d and the scale factor vv
    } // end if
    index[j] = imax;
    // if the pivot is zero the matrix is singular (at least to the precision of
    // the algorithm). For some applications on singular matrices, it is desirable
    // to substitute TINY for zero
    if(0. == a[j][j]) a[j][j] = TINY;
    // now, finally, divide by the pivot element
    if(j != n-1){ // not necessary for the last column
      dum = 1./a[j][j];
      for(int i = j+1; i < n; i++) a[i][j] *= dum;
    } // end if
  } // end for over j, go back for the next column in the reduction
} // end of member function LUDecomposition

/// Solves the equation set ax=b. Here a[0...n-1][0...n-1] is input, not as the
/// matrix A but rather as its LU decomposition, determined by the routine LUDecomposition.
/// index[0...n-1] is input as the permutation vector returned by LUDecomposition.
/// b[0..n-1] is input as the right-hand side and updated with the solution inplace.
/// This routine takes into account the possibility that b will begin with many
/// zeros, so it is efficient for use in matrix inversion
void TAMath::LUBackSubstitution(const matrix &a, int n, const int *index, double *b){
  int ii = -1, ip; double sum;
  for(int i = 0; i < n; i++){
    // exchange b[ip] and b[i] to unscramble the row permutation by LUDecomposition
    sum = b[ip=index[i]];
    if(ip != i) b[ip] = b[i];
    // here begins the forward substitution to solve Ly=b
    if(-1 != ii) for(int j = ii; j < i; j++) sum -= a[i][j]*b[j]; // b[ii] is the 1st non-zero
    else if(sum) ii = i; // element of b, so as to improve calculation efficiency
    b[i] = sum; // confirm the solution of y and stores it in b
  } // end for over i
  for(int i = n; i--;){ // now the backword substitution to solve Ux = y
    sum = b[i];
    for(int j = i+1; j < n; j++) sum -= a[i][j]*b[j];
    b[i] = sum / a[i][i]; // store a component of the slution vector x
  } // end for over i
} // end of member function LUBackSubstitution

/// solve n-dim linar equation set ax=b using LU decomposition
void TAMath::LUSolve(matrix &a, int na, double *b){
  int index[na]; double d;
  LUDecomposition(a, na, index, d);
  LUBackSubstitution(a, na, index, b);
} // end of member function LUSolve
void TAMath::LUSolve(matrix &a, int na, matrix &b, int nb){
  int index[na]; double d;
  LUDecomposition(a, na, index, d);
  // recursively solve b[*][i]
  double bb[na];
  for(int i = nb; i--;){
    for(int j = na; j--;) bb[j] = b[i][j];
    LUBackSubstitution(a, na, index, bb);
    for(int j = na; j--;) b[i][j] = bb[j];
  } // end for over i
} // end of member function LUSolve-1
/// inverse n-dim matrix a using LU decomposition: AX=E
void TAMath::LUInverse(matrix &a, int n, matrix &e){
  int index[n]; double d, b[n];
  if(!a.DimensionMatch(e)) TAException::Error("TAMath", "LUInverse: Dimension Match.");

  LUDecomposition(a, n, index, d);
  for(int j = 0; j < n; j++){ // inverse a column by column - the j-th column
    memset(b, 0, sizeof(b)); b[j] = 1.; // initialize b to the i-th col of the identity matrix
    LUBackSubstitution(a, n, index, b);
    for(int i = n; i--;) e[i][j] = b[i]; // assign b to e[*][j]
  } // end for over i
} // end member function LUInverse
void TAMath::LUInverse(matrix &a, int n){
  matrix e(n, n); // an identity matrix
  for(int i = n; i--;) e[i][i] = 1.;

  LUInverse(a, n, e);
  a = std::move(e);
} // end of member function LUInverse

/// calculate determinent of a square matrix
double TAMath::Det(matrix &a, int n){
  int index[n]; double d;

  LUDecomposition(a, n, index, d);
  for(int i = n; i--;) d *= a[i][i];
  return d;
} // end of member function Det

TEST_CASE("Solve Linear Equation Set", "[linear]"){
  matrix A(3, 3);
  A = { 12., -3.,  3., -18.,  3., -1., 1.,  1.,  1. };
  double b[3] = {15., -15., 6.};
  SECTION("Test Gauss-Jordan Method"){
    TAMath::GaussJordan(A, 3, b);
    CHECK(b[0] == Approx(1.).epsilon(1.e-8));
    CHECK(b[1] == Approx(2.).epsilon(1.e-8));
    CHECK(b[2] == Approx(3.).epsilon(1.e-8));
  }
  SECTION("Test LU Method"){
    TAMath::LUSolve(A, 3, b);
    CHECK(b[0] == Approx(1.).epsilon(1.e-8));
    CHECK(b[1] == Approx(2.).epsilon(1.e-8));
    CHECK(b[2] == Approx(3.).epsilon(1.e-8));
  }
  SECTION("Determinant"){
    CHECK(TAMath::Det(A, 3) == Approx(-66.).epsilon(1.e-8));
  }
} // end of TEST_CASE

/// expand in storage the covariance matrix covar, so as to take into account
/// parameters that are being held fixed, which would be assigned zero covariances
/// \param covar should be a ma*ma 2-D array
/// \param nFit is the count of fitted parameters (not held fixed), -1: not assigned
void TAMath::FillCovar(matrix &covar, int ma, bool *isFit, int nFit){
  if(nFit == -1){ for(int i = nFit = 0; i < ma; i++) if(isFit[i]) nFit++; }
  else if(nFit <= 0) TAException::Error("TAMath", "FillCovar: invalid nFit: %d", nFit);
  for(int i = nFit; i < ma; i++) for(int j = 0; j <= i; j++) covar[i][j] = covar[j][i] = 0.;

  for(int j = ma; j--;){
    if(isFit[j] && j > nFit-1){
      nFit--;
      for(int i = 0; i < ma; i++) swap(covar[i][nFit], covar[i][j]); // col-j <-> col-nFit
      for(int i = 0; i < ma; i++) swap(covar[nFit][i], covar[j][i]); // row-j <-> row-nFit
    } // end if
  } // end for over j
} // end of member function FillCovar
