/**
  \file TADiagonalize.cxx
  \class TADiagonalize
  \brief A collection of matrix diagonalization methods.
  \author SUN Yazhou
  \date Created: 2020/09/27
  \date Last modified: 2020/10/07 by SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cmath>
#include <vector>
#include <algorithm>
#include <catch2/catch.hpp>
#include "TADiagonalize.h"
#include "TAException.h"
#include "TAMatrix.h"
#include "TAMath.h"

using std::vector;
using std::max;
using std::min;

/// computes all eigenvalues and eigenvectors of a real symmetric matrix
/// a[0..n-1][0..n-1]. On output, elements of a are destroyed. d[0..n-1] returns
/// the eigenvalues of a. v is a nxn matrix whose columns contain, on output, the
/// normalized eigenvectors of a. The function returns the number of Jacobi rotations
/// that were required.
/// This implementation is transcribbed from Ref. Numerical Recipes in C, p467.
/// for machines where underflow is not set to zero, the program has to be modified.
inline void rotate(double s, double tau, double &aij, double &akl){
  const double g = aij;
  aij -= s*(akl+aij*tau);
  akl += s*(g  -akl*tau);
} // end rotate
int TADiagonalize::Jacobi(matrix &a, int n, double *d, matrix &v){
  if(!a.IsSymmetric()) TAException::Error("TADiagonalize", "Jacobi: Input matrix not symmetric.");
  int p, q, nrot = 0; // subscripts for the rotation matrix, nrot: the implemented rotations
  const int nn = n, n2 = n*n;
  double b[nn], z[nn]; // b is to accumulate increments(stored in z) of diagonal elements
  double thre, g, h, t; // thre: to choose the povit element in each rotation, g,h,t: temporaries
  // initialize v to identity matrix, and b and d to the diagonal of a //
  for(p = 0; p < n; p++){
    for(q = 0; q < n; q++) v[p][q] = 0.;
    v[p][p] = 1.;
    b[p] = d[p] = a[p][p];
    z[p] = 0.; // accumulate t*a_pq
  } // end for over p
  // here the sweeps begin //
  for(int i = 0; i < 50; i++){
    double sum = 0.; // sum of |off-diagonal| of a
    for(p = 0; p < n-1; p++) for(q = p+1; q < n; q++) sum += fabs(a[p][q]);
    // the normal return, which relies on quadratic convergence to machine precision
    if(0. == sum) return nrot;
    if(i < 3) thre = 0.2*sum/n2; // ...on the first three sweeps
    else thre = 0.; // ...thereafter
    // here begins the sweep - loop over off-diagonal elements //
    for(p = 0; p < n-1; p++){
      for(q = p+1; q < n; q++){
        g = 100.*fabs(a[p][q]);
        // after four sweeps, skip the rotation if the off-diagonal element
        // is too small: less than the least significant digit + 2 of a_pp and a_qq
        if(i > 3 && fabs(d[p]) + g == fabs(d[p]) && fabs(d[q]) + g == fabs(d[q]))
          a[p][q] = 0.;
        else if(fabs(a[p][q]) > thre){
          h = d[q] - d[p];
          if(fabs(h) + g == fabs(h)) t = a[p][q]/h; // t=1/(2theta), for theta is too big
          else{ // calculate t = tan 2phi
            const double theta = 0.5*h/a[p][q];
            t = 1./(fabs(theta)+sqrt(1.+theta*theta));
            if(theta < 0.) t = -t;
          } // end else
          const double c = 1./sqrt(1.+t*t), s = t*c, tau = s/(1.+c);
          h = t*a[p][q];
          z[p] -= h; z[q] += h; // z and d are stored separatly to minimize roundoff
          d[p] -= h; d[q] += h; // d here updated only to help choosing the pivot, i.e. p & q
          a[p][q] = 0.;
          // rotate the above-diagonal elements of a //
          for(int j = 0;   j < p; j++) rotate(s, tau, a[j][p], a[j][q]);
          for(int j = p+1; j < q; j++) rotate(s, tau, a[p][j], a[j][q]);
          for(int j = q+1; j < n; j++) rotate(s, tau, a[p][j], a[q][j]);
          for(int j = 0;   j < n; j++) rotate(s, tau, v[j][p], v[j][q]); // V'=VP
          nrot++;
        } // end else if
      } // end for over q
    } // end for over p
    // update d with the sum of t*apq stored in z, which is stored well in z
    for(p = 0; p < n; p++){
      b[p] += z[p]; // b remained unchanged during the sweep, but now it can accept z[p]
      d[p] = b[p]; // now b contains less roundoff than d, so d is updated with b
      z[p] = 0.;
    } // end for over p
  } // end loop over sweeps
  TAException::Error("TADiagonalize", "Jacobi: Too many sweeps.");
  return nrot; // never gets here
} // end of member function Jacobi

/// note that eigenvalues are not ordered in Jacobi. This routine below sorts
/// the input eigenvalues into descending order, and rearranges the eigenvectors
/// correspondingly. The method is straight insertion
void TADiagonalize::EigenSort(double *d, matrix &v, int n){
  double p;
  int k; // subscript of the maximum
  for(int i = 0; i < n-1; i++){
    p = d[k=i];
    for(int j = i+1; j < n; j++) if(d[j] >= p) p = d[k=j];
    if(k != i){
      d[k] = d[i]; d[i] = p; // exchange di and dk
      for(int j = 0; j < n; j++){ // exchange column i and column k
        p = v[j][i];
        v[j][i] = v[j][k];
        v[j][k] = p;
      } // end for over j
    } // end if
  } // end for over i
} // end of member function EigenSort

/// Householder reduction of a real symmetric matrix to tridiagonal form
/// On output, a is replaced by the orthogonal matrix Q effecting the transformation
/// Note that Q^T=P_1.P_2..P_{n-2}, so it's T=Q^T.a.Q, not the other way around
/// d[0..n-1] returns the diagonal and e[0..n-1] the sub-diagonal, with e[0]=0
/// This implementation is transcribbed from Numerical Recipes in C, p474
#define EIGENVEC // switch on if eigenvectors are wanted
void TADiagonalize::Tridiagonalize(matrix &a, int n, double *d, double *e){
  for(int i = n-1; i >= 1; i--){ // loop over column-i and row-i of a
    int l = i-1; // 0 to l of u are nonzero, the length of u is l+1
    double h = 0., scale = 0.; // scale: for selecting the pivotal column
    if(l > 0){ // nonzeros of u are more than 1
      for(int k = l+1; k--;) scale += fabs(a[i][k]);
      if(0. == scale) e[i] = a[i][l]; // confirm the off-diagonal, and skip the transformation
      else{ // commence the transformation
        for(int k = l+1; k--;){
          a[i][k] /= scale; // scale a to minimize rounding error
          h += TAMath::sqr(a[i][k]); // calculate \sigma (|u|^2) in h
        } // end for over k
        double f = a[i][l]; // the sub-diagonal element of stage i
        double g = (f >= 0. ? -sqrt(h) : sqrt(h)); // g has opposite sign to ai,i-1
        e[i] = scale*g; // +-|x|e1, the off-diagonal, so it is |x| without scale
        h -= f*g; // H=|u|^2/2=x^2-+x1|x|
        a[i][l] = f - g; // store u in the ith row of a: 1st element of u, ai,i-1-+|x|e1
        f = 0.;
        for(int j = 0; j <= l; j++){
#ifdef EIGENVEC // this statement can be omitted if eigenvectors are not wanted
          a[j][i] = a[i][j]/h; // store u/H in ith column of a
#endif
          g = 0.; // form an element of A.u in g
          for(int k =   0; k <= j; k++) g += a[j][k]*a[i][k]; // a is symmetric, rowjA.u
          for(int k = j+1; k <= l; k++) g += a[k][j]*a[i][k]; // only use the sub-diag a, coljA.u
          e[j] = g/h; // p = A.u/H, stored in temporarily unused element of e
          f += e[j]*a[i][j]; // uT.p
        } // end for over j
        const double hh = f/(h+h); // K=uT.p/2H; addition is more efficient than multiplication
        for(int j = 0; j <= l; j++){ // form q and store in e overwriting p
          f = a[i][j]; // uj
          g = e[j] -= hh*f; // q = p - Ku, qj
          for(int k = 0; k <= j; k++) a[j][k] -= f*e[k]+g*a[i][k]; // reduce a, ajk-=uj*qk+uk*qj
        } // end for over j
      } // end else
    } // end if l > 0
    else e[i] = a[i][l]; // no need for reduction for i=1
#ifdef EIGENVEC // this statement can be omitted if eigenvectors are not wanted
    d[i] = h; // =0. if stage-i need not to be reduced (inherently tridiagonal)
#endif
  } // end for over i
  e[0] = 0.;
#ifdef EIGENVEC // this statement can be omitted if eigenvectors are not wanted
  d[0] = 0.; // this statement can be omitted if eigenvectors are not wanted
#endif // this statement can be omitted if eigenvectors are not wanted
  /* Contents of this loop can be omitted if eigenvectors are not wanted
    except for statement d[i] = a[i][i]; */
  for(int i = 0; i < n; i++){
#ifdef EIGENVEC // this statement can be omitted if eigenvectors are not wanted
    if(d[i]){ // i=0 is skipped this way, also skips those stages needing not a transform
      for(int j = 0; j < i; j++){ // d[0] and d[1] are always zero, which insures that...
        double g = 0.; // ...Q_n-2 = P_n-2
        for(int k = i; k--;) g += a[i][k]*a[k][j];
        for(int k = i; k--;) a[k][j] -= g*a[k][i];
      } // end for over j
    } // end if(d[i])
#endif
    d[i] = a[i][i]; // this statement remains
#ifdef EIGENVEC // this statement can be omitted if eigenvectors are not wanted
    // initially Q is E, but since the next iteration only changes the 0~i square patch of Q,
    // and the 0~i-1 square patch has been updated, it suffices to restore the border of the
    // 0~i patch, which are no longer needed, to E. Elements outside the border cannot be
    // altered, because it contains u and u/H that are needed for iterations to follow
    a[i][i] = 1.; // reset row and column of a...
    for(int j = 0; j < i; j++) a[j][i] = a[i][j] = 0.; // ...for next iteration
#endif // this statement can be omitted if eigenvectors are not wanted
  } // end for over i
} // end of member function Tridiagonalize

/// QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors
/// of a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix previously
/// reduced by Tridiagonalize(...) above. On input, d[0..n-1] contains the diagonal
/// elements of the tridiagonal matrix. On output, it returns the eigenvalues. The vector
/// e[0..n-1] inputs the subdiagonal elements of the tridiagonal matrix, with e[0] arbitray.
/// On output e is destroyed. When finding only the eigenvalues, several lines may be omitted,
/// as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired,
/// the matrix z[0..n-1][0..n-1] is input a sthe identity matrix. if the eigenvectors of
/// a matrix that has been reduced by Tridiagonalize are required, then z is input as
/// the matrix output by Tridiagonalize. In either case, the kth column of z returns the
/// normalized eigenvector corresponding to d[k].
void TADiagonalize::TridiagQLImplicit(double *d, double *e, int n, matrix &z){
  for(int i = 1; i < n; i++) e[i-1] = e[i]; // so that e starts from e[0], just for convinence
  e[n-1] = 0.;
  // starts the implicit shifts A-k_l*E //
  // l stands for the order of the l-th eigenvalue, with k_l close to the eigeivalue //
  for(int l = 0; l < n; l++){
    int iter = 0, m;
    do{ // each do-while implements a QL transform, i.e. A'=Q^T.A.Q, until the 2x2
        // block (l,l+1) is diagonalized, (i.e. k_l is used up), then we move to
        // next l, i.e. next k_{l+1}
      // look for a single small subdiagonal element to split the matrix //
      // for matrix that is block-diagonal, to solve the diagonal blocks separately //
      // would give all the eigenvalues and eigenvectors //
      for(m = l; m < n-1; m++){
        const double dd = fabs(d[m])+fabs(d[m+1]);
        if(fabs(e[m])+dd == dd) break;
      } // end for over m
      // here we solve the eigenvalue problems of the the m->l block //
      if(m != l){
        if(iter++ == 30) TAException::Error("TADiagonalize",
          "TridiagQLImplicit: Too many iterations ocurred.");
        // calculate k_l as the eigenvalue of the 2x2 (l,l+1) block //
        double g = (d[l+1]-d[l])/(2.*e[l]);
        double r = TAMath::norm(g, 1.); // sqrt(1.+g*g)
        // g and r have the same sign so that kl is closer to dl, to accelerate convergence
        g = d[m] - (d[l]-e[l]/(g+r*TAMath::sign(g))); // dm-kl
        double s = 1., c = 1., p = 0.; // s,c: the rotation; p: correction to d
        // a plane rotation as in the original QL, followed by Givens rotations
        // to restore tridiagonal form
        int i;
        for(i = m-1; i >= l; i--){ // 1 (m-1,m) elimination + (m-l-1) Givens transform
          // f is the (i-1,i+1) element to be eliminated by current Givens
          double f = s*e[i], b = c*e[i]; // while b is the (i-1,i) element
          e[i+1] = r = TAMath::norm(f, g); // the final e[i+1], afte all the Givens transforms
          if(0. == r){ // find a new split point, exit to start with a new split at m=i+1
            d[i+1] -= p; // complete the unfinished d[i+1] update of the last Givens transform
            e[m] = 0.; break; // recover e[m], and break
          } // end if(0. == r)
          s = f/r; c = g/r; // update s and c for current Givens transform
          g = d[i+1]-p; // update d[i+1] for the last Givens transform in g
          r = (d[i]-g)*s+2.*c*b; d[i+1] = g+(p=s*r); // update d[i+1] for current Givens transform
          g = c*r-b; // updated sub-diagonoal element e[i] after current Givens transform
          // next loop can be omitted if eigenvectors are not wanted //
#ifdef EIGENVEC
          for(int k = 0; k < n; k++){ // accumulate the rotations
            f = z[k][i+1];
            z[k][i+1] = s*z[k][i]+c*f;
            z[k][i] = c*z[k][i]-s*f;
          } // end for over k
#endif
        } // end for over i
        if(0. == r && i >= l) continue; // r=0. represents that the correction to d[i+1] is zero
        d[l] -= p; e[l] = g; //
        e[m] = 0.; // e[m] is spoiled in the loop for convinence, here we recover it
      } // end if(m != l)
    } while (m != l);
  } // end loop over l
} // end of member function TridiagQLImplicit

/// An initial attempt for implementation of Lanczos algorithm ///
/// Lanczos algorithm is devised for a few largest (in modulus) eigenvalues and
/// the corresponding eigenvectors of huge sparse Hermitian matrix A
/// It projects A to a space with much lower dimension (Krylov subspace) as a tridiagonal
/// matrix, whose eigenvalues and eigenvectors are good approximation of A's
/// Given input matrix a[0..n-1][0..n-1], the routine returns first fiew (usually <=10)
/// eigenvalues in d and the corresponding eigenvectors in z
/// Optionally user can provide the initial vector x, or it defaults to {1,1,..,1}
/// d should be of length n, for it is drafted in the program
inline double vprod(const matrix &q, const matrix &r, int n){ // vector inner product
  double p = 0.; for(int i = n; i--;) p += q[i][0]*r[i][0]; // alpha =q*r
  return p;
}
void TADiagonalize::Lanczos(matrix &a, int n, double *d, matrix &z, double *x){
  static const double EPS = 1.e-5; // eigenvalue accuracy standard
  static const int NMAX = 50; // maximum number of iterations
  const int nm = min(n, NMAX), np = nm + 1;
  double e[np]; // d is alpha and e is beta in the conventional notation
  // prepare for the start of the iterations //
  matrix s(nm,nm); // to store the eigenpairs of Tk
  matrix v(n), q(n,1,e), r(n); // q_{j-1}, qj, q_{j+1}, 3 basis vectors
  for(int i = n; i--;) if(x) q[i][0] = x[i]; else q[i][0] = 1.;
  q.cv(0).normalize(); matrix Q(q); // to store the Lanczos vectors
  a.DotProduct(q,r); r.SelfAdd(-(d[0] = vprod(q,r,n)), q); // r = A*q - alphaj*qj
  e[1] = r.cv(0).norm(); // the first off-diagonal of Tk, e[0] is not used
  // OK, here we commence the iteration //
  for(int j = 1; j < nm; j++){
    v = q; // q_{j-1}
    r.Scale(1./e[j],q); // normalize qj in q
    Q.PushBackColumn(q.cv(0));
    a.DotProduct(q, r); // r = A*qj
    r.SelfAdd(-e[j],v); // r -= beta_{j-1}*q_{j-1};
    r.SelfAdd(-(d[j] = vprod(q,r,n)),q); // alpha_j=<qj|A|qj>, r -= alphaj*qj;
    r.Purify(Q); // r -= Q*Q^T*r, to make lanczos vectors more orthogonal XXX
    const int jj = j+1;
    e[jj] = r.cv(0).norm(); // the next beta, equals |r|
    s = 1.; // to store the eigenvector matrix of T_k
    TridiagQLImplicit(d,e,jj,s); // diagonalize T_k
    // calculate the infNorm of betaj*s_{ji}, to determine convergence
    double epsilon = 0.;
    for(int i = jj; i--;) epsilon = max(epsilon, fabs(e[jj]*s[j][i]));
    if(epsilon <= EPS){ Q.DotProduct(s, z); return; } // z=Q*s, the Ritz vector
  } // end for over j
} // end of member function Lanczos

/// Lanczos method featuring Ritz-pairs purging. Unwanted Ritz-pairs are deflated,
/// leaving a smaller Krylov space, where the Lanczos orthogonalization resumes
/// so that the unwanted Ritz-pairs are purged from the final constructed Krylov
/// space. This would speed the convergence to the wanted Ritz-pairs compared with
/// plain Lanczos algorithm. NOTE that we assume the largest eigenvalues are wanted
/// Ref.: http://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf, p212
/// Ref.: K. Wu and H. D. Simon,SIAM J. Matrix Anal. Appl., 22 (2000), pp. 602–616.
void TADiagonalize::LanczosPurge(matrix &a, int n, double *d, matrix &z, double *x){
  static const double EPS = 1.e-8; // eigenvalue accuracy standard
  // number of wanted eigenpairs and dimension for the plain Lanczos
  static const int J = 3, K = J + 2; // J: number of the wanted eigenpairs
  static const int NMAX = 100; // maximum number of iterations
  const int nm = min(n, NMAX);
  const int M = max(nm, K+4), np = nm+1; // dimension of the final Tk
  double e[np]; // d is alpha and e is beta in the conventional notation
  // prepare for the start of the iterations //
  matrix v(n), q(n), r(n); // q_{j-1}, qj, q_{j+1}, 3 basis vectors
  matrix s(M,M), st(M,M); st = 1.; // s: the eigenvector matrix of Tk, st set to E
  for(int i = n; i--;) if(x) q[i][0] = x[i]; else q[i][0] = 1.;
  q.cv(0).normalize(); matrix Q(q); // to store the Lanczos vectors
  a.DotProduct(q,r); r.SelfAdd(-(d[0] = vprod(q,r,n)), q); // r = A*q - alphaj*qj
  e[1] = r.cv(0).norm(); // the first off-diagonal of Tk, e[0] is not used
  // commence the iteration
  for(int j = 1; j < nm; j++){
    v = q; // q_{j-1}
    r.Scale(1./e[j],q); // normalize qj in q
    Q.PushBackColumn(q.cv(0));
    a.DotProduct(q,r); // r = A*qj
    r.SelfAdd(-e[j],v); // r -= beta_{j-1}*q_{j-1};
    r.SelfAdd(-(d[j] = vprod(q,r,n)),q); // alpha_j=<qj|A|qj>, r -= alphaj*qj;
    r.Purify(Q); // r -= Q*Q^T*r, to make lanczos vectors more orthogonal XXX
    const int jj = j+1;
    e[jj] = r.cv(0).norm(); // betaj=|r|
    if(K+1 == j) TridiagArrow(d,e,J+1,st); // tridiagonalize the mutant Tk
    s = st; // acknowledge the transform from arrow to tridiagonal
    // diagonalize T_k //
    TridiagQLImplicit(d,e,jj,s); // eigenvalues in descending order is mandatory
    // calculate the infNorm of betaj*s_{ji}
    double epsilon = 0.;
    for(int i = jj; i--;) epsilon = max(epsilon, fabs(e[jj]*s[j][i]));
    if(epsilon <= EPS){ Q.DotProduct(s,z); return; } // z=Q*s, the Ritz vector
    if(K-1 == j){ // K Ritz-vecs obtained, initiate the restart
      // keep only the wanted Ritz-pairs, and update the tridiagonal matrix //
      r.Scale(1./e[jj],q); // q_{k+1} = r/|r|
      a.DotProduct(q,r); // r = a*q_{k+1}
      r.SelfAdd(-d[j],q); // remove the q_{k+1} component
      Q.DotProduct(s,z); // z=Q*s: calculate the Ritz vector
      Q.FreeMemory(); // make space for the Ritz vectors
      for(int i = 0; i < J; i++){
        Q.PushBackColumn(z.cv(i)); // confirm zi as the new basis, update Tk and..
        r.SelfAdd(-(e[i+1] = e[jj]*s[j][i]), z.cv(i)); // ..remove zi bits from r
      } // end for over i
      // OK, here we have the starting vector now, and the plain Lanczos..
      d[jj] = r.cv(0).norm(); // ..can be restarted from here
      j -= K-J; n -= K-J; // one more thing, squeeze the loop for the purge
    } // end if(j == K)
  } // end for over j
} // end of member function LanczosRestart

TEST_CASE("Test diagonalization algorithm", "[diag]"){
  static const int n = 5;
  matrix a(n,n);
  a = {1,2,3,4,5, 2,3,4,5,1, 3,4,5,1,2, 4,5,1,2,3, 5,1,2,3,4};
  double d[n];
  SECTION("Jacobi"){
    TADiagonalize::JacobiSort(a,n,d);
    // matrix(1,n,d).Print(); a.Print();
    CHECK(d[0] == Approx(15.).epsilon(1.e-10));
    CHECK(d[1] == Approx(4.25325404176).epsilon(1.e-10));
    CHECK(d[2] == Approx(2.62865556059566).epsilon(1.e-10));
    CHECK(a[0][0] == Approx(0.4472135955).epsilon(1.e-10));
    CHECK(a[3][2] == Approx(0.563522005301).epsilon(1.e-10));
    CHECK(a[4][3] == Approx(-0.563522005301).epsilon(1.e-10));
  } // end SECTION Jacobi
  SECTION("Tridiagonalize"){
    TADiagonalize::TridiagQLSort(a,n,d); // more efficient than JacobiSort
    CHECK(d[0] == Approx(15.).epsilon(1.e-10));
    CHECK(d[1] == Approx(4.25325404176).epsilon(1.e-10));
    CHECK(d[2] == Approx(2.62865556059566).epsilon(1.e-10));
    CHECK(a[0][0] == Approx(0.4472135955).epsilon(1.e-10));
    CHECK(a[3][2] == Approx(-0.563522005301).epsilon(1.e-10));
    CHECK(a[4][3] == Approx(-0.563522005301).epsilon(1.e-10));
  } // end SECTION Tridiagonalize
  SECTION("multiplicity"){
    double dd[] = {1,2,2,3,4,5}; a.diag(dd);
    TADiagonalize::TridiagQLSort(a,n,d);
    CHECK(d[0] == Approx(4.).epsilon(1.e-10));
    CHECK(d[2] == Approx(2.).epsilon(1.e-10));
    CHECK(d[3] == Approx(2.).epsilon(1.e-10));
  } // end SECTION multiplicity
} // end of TEST_CASE

TEST_CASE("Lanczos Test", "[lancz]"){
  static const int n = 6;
  matrix a(n,n), z(n,n); // automatically all initialized to 0
  double d[n] = {0,1,2,3,4,100000}, x[n] = {1,1,1,1,1,1}; // (z,d) are the eigenpair
  for(int i = n; i--;) a[i][i] = d[i];
  a.Print();
  TADiagonalize::Lanczos(a,n,d,z,x);
} // end of TEST_CASE lanczos
