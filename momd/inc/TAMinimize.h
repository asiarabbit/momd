/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMinimize.h
  \class TAMinimize
  \brief toolkit for multidimensional function minimizations. This is supposed to
  be a tool, thus static class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/09/16
  \date Last modified: 2020/09/17 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMinimize_h_
#define _TAMinimize_h_

template <typename FUNC>
class TAMinimize{
public:
  TAMinimize(){}
  virtual ~TAMinimize(){}

  /// to bracket a minimum point of function f in an interval (a,b,c), where
  /// f(a)>f(b) and f(c)>f(b). a,b,c and their respective function values fa, fb
  /// and fc are stored. a and b have to be given upon calling this method
  /// FUNC has to be double (*f)(double) or has member function double operator(double)
  static void Bracket(double &a, double &b, double &c,
    double &fa, double &fb, double &fc, FUNC &f);

  /// given a function f, and a bracketing triplet of abscissas a, b, c such that
  /// b is between a and c, and f(b) is less than both f(a) and f(b), this method
  /// performs a golden section search for the minimum, isolating it to a fractional
  /// precision of about tol. The abscissas of the minimum is returned as xmin, and
  /// the minimum is returned as the return value of the method
  /// note that tol should be larger than sqrt(machine precision of double)~3e-8,
  /// or the function value difference is smaller than roundoff error.
  /// FUNC has to be double (*f)(double) or has member function double operator(double)
  static double Golden(double a, double b, double c, FUNC &f, double &xmin,
      double tol = 1.e-6);

  /// Brent's minimization scheme, which employs golden section and parabolic interpolation
  /// in a coorperative manner. Given a function f, a bracketing triplet abscissas ax, bx, and cx
  /// such that bx is between ax and cx, and f(x) is less than both f(ax) and f(cx), this routine
  /// isolates the minimum to a fractional precision of about tol. The minimum is returned as the
  /// returning value, and the corresponding abscissa is stored in xmin
  /// note that tol should be larger than sqrt(machine precision of double)~3e-8,
  /// or the function value difference is smaller than roundoff error.
  /// FUNC has to be double (*f)(double) or has member function double operator(double)
  static double Brent(double ax, double bx, double cx, FUNC &f, double &xmin,
      double tol = 1.e-6);
  static double Brent(double a, double b, FUNC &f, double &xmin, double tol = 1.e-6){
    double c, fa, fb, fc;
    Bracket(a, b, c, fa, fb, fc, f); // determine a triplet (a,b,c) to wrap a minimum in
    return Brent(a,b,c,f,xmin,tol);
  }
  /// the version of brent using the function's first derivative. df is the first derivative
  static double DBrent(double ax, double bx, double cx, double (*f)(double), double (*df)(double),
      double &xmin, double tol = 1.e-6);

  /// n-dimensional minimization routine using downhill simplex method
  /// \param p stores (n+1) vectors of dimension n. Upon returning, |fmax-fmin|/(|fmax|+|fmin|)<2tol
  /// \param nf tracks the number of function evaluation (calls of f)
  /// NOTE that y must be preinitialized to f(p)
  /// FUNC has to be double (*f)(double *) or has member function double operator(double*)
  static void Amoeba(matrix &p, double *y, int n, double ftol, FUNC &f, int &nf);
  /// extrapolate by a factor fac through the face of the simplex across from the high point, tries it, and
  /// replaces the high point if the new point is better
  static double Amotry(matrix &p, double *y, double *psum, int n, FUNC &f, int ih, double fac);
};

#include "TAMinimize.hpp"

#endif
