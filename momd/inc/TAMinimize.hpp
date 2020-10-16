/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMinimize.hpp
  \class TAMinimize
  \brief toolkit for multidimensional function minimizations.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/09/17
  \date Last modified: 2020/10/05 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <algorithm>
#include "TAException.h"
#include "TAMath.h"

#define sign TAMath::sign

using std::max;
using std::swap;

/// to bracket a minimum point of function f in an interval (a,b,c), where
/// f(a)>f(b) and f(c)>f(b). a,b,c and their respective function values fa, fb
/// and fc are stored. a and b have to be given upon calling this method
static const double GOLD = TAMath::Alpha(); // 0.61803398874989484820458683436565
static const double GOLDP = 1. + GOLD; // 1.61803398874989484820458683436565
static const double GOLDM = 1. - GOLD; // 1-0.61803398874989484820458683436565
static const double TINY = 1.e-20; // to prevent any possible division by zero
static const double GLIM = 100.; // to prevent any possible division by zero
inline void shift(double &a, double &b, double &c, double d){ a = b; b = c; c = d; }
inline void shift(double &a, double &b, double c){ a = b; b = c; }
inline double gold(double &a, double &b){ return b + GOLDP*(b-a); } // inrement the way from a to b
template <typename FUNC>
void TAMinimize<FUNC>::Bracket(double &a, double &b, double &c,
    double &fa, double &fb, double &fc, FUNC &f){ // f: double (*f)(double)
  fa = f(a); fb = f(b); // evaluate f(a) and f(b)
  if(fa < fb){ swap(a, b); swap(fa, fb); } // make sure that fa>=fb, so that a to be is downhill
  // expand original interval (a,b) for an interval (a,b,c) where f(a)>f(b) and f(c)>f(b)
  // Why golden section ratio?
  // Because the incremented interval is equal to the original search domain
  // i.e. for (a,b,c) and (b,c,c'), |c'-c|=|c-a|, so as to accelerate the searching process
  // and at the same time, |b-c|/|b-a|=|c'-c|/|c-b|, i.e. scale similarity achieved
  c = gold(a, b); fc = f(c); // first guess for c
  while(fb > fc){ // keep returning here until we bracket
    // fit function f using (a,fa), (b,fb) and (c,fc) to compute the minimum point u of the parabola
    const double r = (b-a)/(fb-fc);
    const double q = (b-c)/(fb-fa);
    double u = b - (q*(b-c)-r*(b-a)) /
      (2.*max(fabs(q-r), TINY)*sign(q-r)), fu; // to prevent any possible division by zero
    const double ulim = b+GLIM*(c-b); // the limit for u, we won't go farther than this
    // test various possibilities //
    if((b-u)*(u-c) > 0.){ // parabolic u is between b and c: try it
      fu = f(u);
      if(fu < fc){ // fb>fc>fu, i.e. fb>fu and fc>fu. There's got to be a minimum in (b,u,c)
        // so we return (b,u,c) to (a,b,c) //
        a  = b; fa = fb;
        b  = u; fb = fu;
        return;
      } // end if(fu < fc)
      if(fu > fb){ // fa>fb, fu>fb, so a minimum in (a,b,u) is guaranteed.
        c = u; fc = fu; return; // we return (a,b,u) to (a,b,c)
      } // end if(fu > fb)
      // parabolic fit was no use. Use default magnification and start over in while. Very rare
      u = gold(b, c); fu = f(u); // then (b,c,u)->(a,b,c), in the end of while
    } // end if u is between b and c
    if((c-u)*(u-ulim) > 0.){ // u is between c and ulim
      fu = f(u);
      if(fu < fc){ // still downhill, so the 3rd point is beyond u,
        // advance a further step, then (b,c,u)->(a,b,c), in the end of while
        shift(b,  c,  u,  gold(c, u));
        shift(fb, fc, fu, f(u));
      } // or fu>=fc, then (b,c,u) brackets, then (b,c,u)->(a,b,c), in the end of while
    }
    if((u-ulim)*(ulim-c) >= 0.){ // u exceeds ulim,
      u = ulim; fu = f(u); // then we pull it back to ulim, and keep searching
    }
    else{ // reject parabolic u, and use default magnification. Very rare
      u = gold(b, c); fu = f(u); // then (b,c,u)->(a,b,c), in the end of while
    }
    // do the assignment (b,c,u)->(a,b,c), and start the search over again
    shift(a,b,c,u);
    shift(fa,fb,fc,fu);
  } // end while
} // end of member function Bracket

/// given a function f, and a bracketing triplet of abscissas a, b, c such that
/// b is between a and c, and f(b) is less than both f(a) and f(b), this method
/// performs a golden section search for the minimum, isolating it to a fractional
/// precision of about tol. The abscissas of the minimum is returned as xmin, and
/// the minimum is returned as the return value of the method
/// note that tol should be larger than sqrt(machine precision of double)~3e-8,
/// or the function value difference is smaller than roundoff error.
template <typename FUNC>
double TAMinimize<FUNC>::Golden(double a, double b, double c, FUNC &f, // double (*f)(double)
    double &xmin, double tol){
  if((a-b)*(b-c) < 0.) TAException::Error("TAMinimize",
    "Golden: b(%f) is not between a(%f) and c(%f).", b, a, c);
  double x0 = a, x1, x2, x3 = c; // the four points to contract the bracketing interval
  // insert x2 into the larger segment //
  if(fabs(c-b) > fabs(b-a)){ x1 = b; x2 = x3 + GOLD*(x1-x3); } // insert x2 in (x1,x3)
  else{ x2 = b; x1 = x0 + GOLD*(x2-x0); } // |b-c| is smaller, then assign x2 to b,
  double f1 = f(x1), f2 = f(x2); // and
  // start the searching loop. f is called only once for each loop
  while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))){
    // we choose the new trial point so that the minimum one is never at either of the ends
    if(f2 < f1){ // wrap x2 in
      shift(x0, x1, x2, x3+GOLD*(x2-x3)); // x2-x3: the new x1-x3
      shift(f1, f2, f(x2));
    } // end if
    else{
      shift(x3, x2, x1, x0+GOLD*(x1-x0)); // x1-x0: the new x2-x0
      shift(f2, f1, f(x1));
    } // end else
  } // end while
  // now with the desired accuracy, we return the minimum //
  if(f1 < f2){ xmin = x1; return f1; }
  else{ xmin = x2; return f2; }
} // end of member function Golden

/// Brent's minimization scheme, which employs golden section and parabolic interpolation
/// in a coorperative manner. Given a function f, a bracketing triplet abscissas ax, bx, and cx
/// such that bx is between ax and cx, and f(x) is less than both f(ax) and f(cx), this routine
/// isolates the minimum to a fractional precision of about tol. The minimum is returned as the
/// returning value, and the corresponding abscissa is stored in xmin
/// note that tol should be larger than sqrt(machine precision of double)~3e-8,
/// or the function value difference is smaller than roundoff error.
template <typename FUNC> // double (*f)(double)
double TAMinimize<FUNC>::Brent(double ax, double bx, double cx, FUNC &f, double &xmin,
    double tol){
  if((ax-bx)*(bx-cx) < 0.) TAException::Error("TAMinimize",
    "Brent: bx(%f) is not between ax(%f) and cx(%f).", bx, ax, cx);
  static const int ITMAX = 100; // maximum iteration times
  static const double ZEPS = 1.e-10; // protects against division by 0
  double a, b, x, w, v, u; // w,v, and u in (a,x,b), where (fx,fw,fv) in ascending order
  double fx, fw, fv, fu; // function value of x,w,v,u, repsectively
  a = ax < cx ? ax : cx; // so that a > b
  b = ax < cx ? cx : ax;
  x = w = v = bx; fx = fw = fv = f(v); // initialization
  double etmp, e = 0., d; // e: the step before last, d: the last step
  for(int i = 0; i < ITMAX; i++){
    const double xm = 0.5*(a+b); // the mid point of (a,b)
    const double tol1 = tol*fabs(x)+ZEPS, tol2 = 2.*tol1;
    if(0.5*(b-a)+fabs(x-xm) <= tol2){ xmin = x; return fx; } // |x-xm|<|b-a|/2<tol*|x|=>|x-xm|/x<tol
    // compute new point u: golden section or parabolic interpolation
    if(fabs(e) > tol1){
      double q = (x-v)*(fx-fw), r = (x-w)*(fx-fv);
      double p = (x-v)*q - (fx-fw)*r;
      q = 2.*(q-r);
      if(q > 0.) p = -p; // x' = x0+p/q = x0+d, where d = p/q
      q = fabs(q);
      shift(etmp, e, d); // save the previous e and d
      if(fabs(p) > fabs(0.5*q*etmp) || p <= q*(a-x) || p >= q*(b-x)) // d > 0.5etmp or x outside (a,b)
        d = GOLDM*(e = x >= xm ? a-x : b-x); // golden section of (a,x) or (b,x), the longer one
      else{ // use parabolic interpolation
        u = x + (d = p/d); // the parabolic step
        if(u-a < tol2 || b-u < tol2) // if u is < 2*tol from the either end of (a,b),
          d = sign(d)*tol1; // original code in NR(p405): TAMath::sign(xm-x)*tol1: FIXME
      } // end else
    } // end if(|e|>tol1)
    else d = GOLDM*(e = x >= xm ? a-x : b-x); // golden section of (a,x) or (b,x), the longer one
    // OK, now we get the step //
    u = fabs(d) >= tol1 ? x + d : x + sign(d)*tol1; fu = f(u); // the only call of f per iteration
    // now we arrange to ensure that (fx,fw,fv) are in ascending order
    if(fu <= fx){
      if(u >= x) a = x; else b = x; // make x the new border, and wrap u in
      shift(v,w,x,u); shift(fv,fw,fx,fu);
    } // end if(fu >= fx)
    else{
      if(u < x) a = u; else b = u; // make u the new border, and wrap x in
      // fx<=fu<=fw, or w==x (fw=fx<=fu), u is the 2nd minimum, assign w to it
      if(fu <= fw || w == x){ shift(v,w,u); shift(fv,fw,fu); }
      else if(fu <= fv || v == x || v == w){ v = u; fv = fu; } // fx<fw<fu, fu the 3rd min, let v->u
      // the last scenario, if(fu > fw and fu>fv), do nothing. The parabolic triplet remains //
    } // end else
  } // end for over i
  TAException::Error("TAMinimize", "Brent: Too many iterations occurred.");
  xmin = x; return fx; // never gets here
} // end of member function Brent

/// the version of brent using the function's first derivative. df is the 1st derivative
template <typename T>
inline void move(T &a, T &b, T &c,   const T &d, const T &e, const T &f){ a = d; b = e; c = f; }
template <typename FUNC>
double TAMinimize<FUNC>::DBrent(double ax, double bx, double cx, double (*f)(double), double (*df)(double),
    double &xmin, double tol){
  if((ax-bx)*(bx-cx) < 0.) TAException::Error("TAMinimize",
    "Brent: bx(%f) is not between ax(%f) and cx(%f).", bx, ax, cx);
  static const int ITMAX = 100; // maximum iteration times
  static const double ZEPS = 1.e-10; // protects against division by 0
  double a, b, d, u, v, w, x, d1, d2, e = 0.;
  double fu, fv, fw, fx, du, dv, dw, dx; // f & f' of x,w,v,u, repsectively
  a = ax < cx ? ax : cx; // so that a > b
  b = ax < cx ? cx : ax;
  x  = w  = v  = bx;
  fx = fw = fv = f(v);
  dw = dv = dx = df(x); // initialization
  for(int i = 0; i < ITMAX; i++){
     const double xm = 0.5*(a+b);
     const double tol1 = tol*fabs(x)+ZEPS, tol2 = 2.*tol1;
     if(0.5*(b-a)+fabs(x-xm) <= tol2){ xmin = x; return fx; }
     if(fabs(e) > tol1){
       double d1 = 2.*(b-a), d2 = d1;
       if(dw != dx) d1 = (w-x)*dx/(dx-dw); // secant method with one point
       if(dv != dx) d2 = (v-x)*dx/(dx-dv); // and the other
       const double u1 = x + d1, u2 = x + d2;
       const bool ok1 = (a-u1)*(u1-b) > 0. && dx*d1 <= 0.; // u within (a,b) and
       const bool ok2 = (a-u2)*(u2-b) > 0. && dx*d2 <= 0.; // dx*d1(2) <= 0 so that d1(2) is downhill
       const double olde = e; // movement on the step before last
       e = d;
       if(ok1 || ok2){ // take only an acceptable d, and if both are acceptable,
         if(ok1 && ok2) d = fabs(d1) < fabs(d2) ? d1 : d2; // then take the smallest one
         else if(ok1) d = d1; else d = d2;
         if(fabs(d) <= fabs(0.5*olde)){ // current step is half the step before last
           u = x + d;
           if(u-a < tol2 || b-u < tol2) d = sign(d)*tol1; // falls in the forbidden region, truncate d
         } // end if(|d|<|0.5*olde|)
         else d = 0.5*(e = dx >= 0. ? a-x : b-x); // decide which segment by the sign of the derivative
       } // end if(ok1||ok2)
       else d = 0.5*(e = dx >= 0. ? a-x : b-x); // bisect, not golden section (rootfinding of f')
     } // end if(|e| > tol1)
     else d = 0.5*(e = dx >= 0. ? a-x : b-x);
     // here we've got the next step d
     if(fabs(d) >= tol1){ u = x + d; fu = f(u); }
     else{
       u = x + sign(d)*tol1;
       fu = f(u);
       // if the minimum step in the downhill direction take us uphill, then we are done
       if(fu > fx){ xmin = x; return fx; }
     } // end else
     // proceed with step d //
     du = df(u);
     // prepare for the next iteration. Everything is the same as method Brent //
     if(fu <= fx){ // push x to the border, and then assign u to it
       if(x <= u) a = x; else b = x;
       move(v, fv, dv,  w, fw, dw);
       move(w, fw, dw,  x, fx, dx);
       move(x, fx, dx,  u, fu ,du);
     } // end if(fu <= fx)
     else{ // push u to the border
       if(u < x) a = u; else b = u;
       if(fu < fw || w == x){ // (fx,fu,fw), assign u to w
         move(v, fv, dv,  w, fw, dw);
         move(w, fw, dw,  u, fu, du);
       } // end if fu is the 2nd largest
       else if(fu < fv || v == x || v == w){ // (fx,fw,fu,fv), fu is the 3rd largest, v->u
         move(v, fv, dv,  u, fu, du);
       } // end elseif
     } // end else
  } // end for over i
  TAException::Error("TAMinimize", "DBrent: Too many iterations occurred.");
  return 0.; // never get here
} // end of member function DBrent


/// n-dimensional minimization routine using downhill simplex method
/// \param p stores (n+1) vectors of dimension n. Upon returning, |fmax-fmin|/(|fmax|+|fmin|)<2tol
/// \param nf tracks the number of function evaluation (calls of f)
/// NOTE that y must be preinitialized to f(p)
template <typename FUNC> // double (*f)(double *)
void TAMinimize<FUNC>::Amoeba(matrix &p, double *y, int n, double ftol, FUNC &f, int &nf){
  static const double TINY = 1.e-10; // a small number
  static const double NMAX = 5000; // a small number
  if(p.nr() != n+1 || p.nc() != n) TAException::Error("TAMinimize", "Amoeba: Dimension mismatch.");
  nf = 0;
  const int nn = n, np = n+1;
  double psum[nn], sum;
  int i, j, ih, inh, il; // the point of the highest, next-highest, and the lowest
  for(j = n; j--;){
    for(i = np, sum = 0.; i--;) sum += p[i][j];
    psum[j] = sum;
  } // end for over j
  while(1){
    //// TEST FOR CONVERGENCE FIRST ////
    il = 0; // the point with lowest function value
    // first we must determine which point is the highest(worst, next-highest,
    // and lowest(best), by looping over the points in the simplex
    int ih = y[0] > y[1] ? (inh=1, 0) : (inh=0, 1);
    for(i = np; i--;){
      if(y[i] <= y[il]) il = i;
      if(y[i] > y[ih]){ inh = ih; ih = i; }
      else if(y[i] > y[inh] && i != ih) inh = i;
    } // end for over i
    // compute the fractional range from highest to lowest and return if satisfactory
    const double rtol = 2.*fabs(y[ih]-y[il])/(fabs(y[ih])+fabs(y[il])+TINY);
    if(rtol < ftol){ // if returning, put the best point and value in slot 0
      swap(y[0], y[ih]);
      for(i = n; i--;) swap(p[0][i], p[il][i]);
      break;
    } // end if(rtol<ftol)
    if(nf >= NMAX) TAException::Error("TAMinimize", "Amoeba: Too many iterations occurred.");
    // start going simplex-downhill //
    nf += 2;
    // first, reflect: extrapolate by a factor -1 through the face of the simplex
    // across from the high point, i.e., reflect the simplex from the high point
    double ytry = Amotry(p, y, psum, n, f, ih, -1.);
    if(ytry <= y[il]) ytry = Amotry(p, y, psum, n, f, ih, 2.); // Good, advance another same step
    else if(ytry >= y[inh]){ // worst than the second-highest, let's retract by half a step
      const double ysave = y[ih];
      ytry = Amotry(p, y, psum, n, f, ih, 0.5);
      if(ytry >= ysave){ // still no signs of downhill. Better contract around point_il
        for(i = np; i--;) if(i != il){
          for(j = n; j--;) p[i][j] = psum[j] = 0.5*(p[il][j]+p[i][j]);
          y[i] = f(psum); // totally n calls of f
        } // end for over i
      } // end if(ytry >= ysave)
      nf += n; // keep track of function evaluations
      for(j = n; j--;){ // update psum
        for(i = np, sum = 0.; i--;) sum += p[i][j];
        psum[j] = sum;
      } // end for over j - end update psum
    } // end if(ytry >= y[inh])
    else nf--; // correct the evaluation count
  } // end while
} // end of member function Ameoba
/// extrapolate by a factor fac through the face of the simplex across from the high point, tries it, and
/// replaces the high point if the new point is better
template <typename FUNC> // double (*f)(double *)
double TAMinimize<FUNC>::Amotry(matrix &p, double *y, double *psum, int n, FUNC &f,
    int ih, double fac){
  const int nn =  n; double ptry[nn], ytry;
  const double fac1 = (1.-fac) / n, fac2 = fac1 - fac;
  for(int j = n; j--;) ptry[j] = psum[j]*fac1 - p[ih][j]*fac2; // the reflected point
  ytry = f(ptry); // evaluate the function at the trial point
  if(ytry < y[ih]){ // accept ptry, replace the highest
    y[ih] = ytry;
    for(int j = n; j--;){
      psum[j] += ptry[j] - p[ih][j];
      p[ih][j] = ptry[j];
    }
  } // end if
  return ytry;
} // end of member function Amotry
