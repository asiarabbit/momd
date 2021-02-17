/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TALanczos.cxx
  \class TALanczos
  \brief Lanczos algorithm for solving large sparse matrix iteratively. Transcribbed
  from momd/TADiagonalize::Lanczos and momd/TADiagonalize::LanczosPurge, but with
  new classes representing large sparse matrix using ROOT - TTree object storing
  large sparse matrix and vectors, so as to save memory. This is a static class.
  Dedicated for diagonalization problems arising from large-scale shell model calculations.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/12, Chinese Lunar New Year
  \date Last modified: 2021/02/12 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <vector>
#include <algorithm>
#include "TALanczos.h"
#include "TADiagonalize.h"
#include "TAException.h"
#include "TAHamiltonian.h"
#include "TASparseVec.h"
#include "TASparseMatrix.h"
#include "TRandom3.h"

using std::vector;
using std::min;
using std::max;

typedef unsigned long long ull;

/// Lanczos algorithm is devised for a few largest (in modulus) eigenvalues and
/// It projects A to a space with much lower dimension (Krylov subspace) as a tridiagonal
/// matrix, whose eigenvalues and eigenvectors are good approximation of A's
/// Given input matrix a[0..n-1][0..n-1], the routine returns first fiew (usually <=10)
/// eigenvalues in d and the corresponding eigenvectors in z
/// d should be at least of length 50, for it is drafted in the program
void TALanczos::Lanczos(TAHamiltonian &H, double *d, TASparseMatrix &z, TASparseVec *x, int nPair){
  static const double EPS = 1.e-6; // eigenvalue accuracy standard
  static const int NMAX = 50; // maximum number of iterations
  const int nm = max(nPair, NMAX), np = nm + 1;
  double e[np]; vector<double> alpha, beta;
  // prepare for the start of the iterations //
  TASparseVec v("v"), q("q"), r("r"); // q_{j-1}, qj, q_{j+1}, 3 basis vectors
  if(x && x->GetEntries()) q = *x; else q.Fill(0, 1.); // initialize q to (1,0,0,..)
  q.normalize();
  TASparseMatrix Q; Q.PushBackColumn(q); // to store the Lanczos vectors
  H.DotProduct(q, r); // r = H*q
  alpha.push_back(r.DotProduct(q)); // alpha = <qj|a|qj>
  r.SelfAdd(-alpha[0], q); // r -= alphaj*qj
  beta.push_back(r.norm()); // the first off-diagonal of Tk
  // the trial vector is an eigenvector
  if(fabs(beta[0]) < EPS){
    d[0] = alpha[0];
    z.Clear(); z.PushBackColumn(q);
    z.Print();
    return;
  }
  // OK, here we commence the iteration //
  for(int j = 1; j <= nm; j++){
    v = q; // q_{j-1}
    r.Scale(1./beta[j-1],q); // normalize qj in q
    Q.PushBackColumn(q);
    H.DotProduct(q, r); // r = H*qj
    r.SelfAdd(-beta[j-1], v); // r -= beta_{j-1}*q_{j-1};
    alpha.push_back(r.DotProduct(q)); // alpha_j=<qj|a|qj>
    r.SelfAdd(-alpha[j], q); // r -= alphaj*qj;
    r.Purify(Q); // r -= Q*Q^T*r, correct roundoff error for loss of orthogonality
    beta.push_back(r.norm()); // the next beta, equals |r|
    const int jj = j+1;
    for(int i = jj; i--;) d[i] = alpha[i]; // diagonal of Tk
    for(int i = j;  i--;) e[i+1] = beta[i]; // e[0] is not used
    matrix s(jj,jj); s = 1.; // to store the eigenvector matrix of Tk
    TADiagonalize::TridiagQLImplicit(d,e,jj,s); // diagonalize Tk
    double epsilon = 0.; // check convergence
    for(int i = jj; i--;) epsilon = max(epsilon, fabs(beta[j]*s[j][i]));
    if(epsilon <= EPS || nm == j){
      TADiagonalize::EigenSort(d,s,jj);
      Q.DotProduct(s, z);
      // r.Print(); r.Print(); r.Save(); return; // DEBUG
      if(nm == j) break; else return;
    }
  } // end for over j
  TAException::Info("TALanczos", "Lanczos: Too many iterations occurred.");
} // end of member function Lanczos

/// Lanczos method featuring Ritz-pairs purging. After getting a set of Ritz vectors
/// from a plain Lanczos, we keep a smaller set and remove (purge) the rest. Then
/// develop Krylov space from the new Lanczos basis. So unwanted directions are removed
/// from the initial vector for each iteration. NOTE that the largest eigenvalues are kept.
/// \param nPair: number of eigenpairs wanted. Still the function will return as long
/// \param x: the initial vector to initiate the iteration
/// as encountering an invariant space, even before finding nPair eigenpairs. For this
/// drawback, users are recommended with member method PurgeRestart(...)
/// Ref.: http://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf, p212
/// Ref.: K. Wu and H. D. Simon,SIAM J. Matrix Anal. Appl., 22 (2000), pp. 602–616.
void TALanczos::LanczosPurge(TAHamiltonian &H, double *d, TASparseMatrix &z, TASparseVec *x, int nPair){
  static const double EPS = 1.e-8; // eigenvalue accuracy standard
  static const int NMAXRST = 100; // maximum restart times
  // number of wanted eigenpairs and dimension for the plain Lanczos
  const int J = nPair, K = J + 2; // J, K: number of the wanted eigenpairs and candidates
  const ull n = H.GetNBasis();
  if(n <= ull(K)) return Lanczos(H,d,z,x,nPair); // we deal with a little bit larger matrix
  const int NMAX = 20; // maximum number of iterations
  const int nm = min(n, (ull)NMAX);
  const int np = nm+1; // maximum dimension of the final Tk
  double e[np]{}; vector<double> alpha, beta;
  int restartCnt = 0; // restart times
  bool restarted = false; // if has just restarted
  // prepare for the start of the iterations //
  matrix st(J+1,J+1); st = 1.; // s: the eigenvector matrix of Tk, st set to E
  TASparseVec v("v"), q("q"), r("r"); // q_{j-1}, qj, q_{j+1}, 3 basis vectors
  if(x && x->GetEntries()) q = *x; else q.Fill(0, 1.); // initialize q to (1,0,0,..)
  q.normalize();
  TASparseMatrix Q; Q.PushBackColumn(q); // to store the Lanczos vectors
  H.DotProduct(q, r); // r = a*q
  alpha.push_back(r.DotProduct(q)); // alpha = <qj|a|qj>
  r.SelfAdd(-alpha[0], q); // r = A*q-alphaj*qj
  beta.push_back(r.norm()); // the first off-diagonal of Tk
  // FIXME: if the trial vector is an eigenvector, return
  if(fabs(beta[0]) < EPS){
    d[0] = alpha[0];
    z.Clear(); z.PushBackColumn(q);
    return;
  }
  // commence the iteration //
  for(int j = 1; j <= nm; j++){ // j: n of Lanczos vectors, so it'd be always far from nm
    v = q; // q_{j-1}
    r.Scale(1./beta[j-1], q); // normalize r in q
    H.DotProduct(q, r); // r = a*qj
    r.SelfAdd(-beta[j-1], v); // r -= beta_{j-1}*q_{j-1};
    alpha.push_back(r.DotProduct(q)); // alpha=<qj|a|qj>
    r.SelfAdd(-alpha[j], q); // r -= alphaj*qj;
    const int jj = j+1;
    // restoring d, e and s is nontrivial, for they are changed per iteration
    for(int i = jj; i--;) d[i] = alpha[i]; // diagonal of Tk
    for(int i = j; i--;) e[i+1] = beta[i]; // e[0] is not used
    if(restarted){ // just restarted
      TADiagonalize::TridiagArrow(d,e,j,st); // tridiagonalize the mutant Tk to extract the Lanczos basis
      for(int i = j; i--;) alpha[i] = d[i]; // Tk recovered to tridiagonal, j=J+1
      for(int i = J; i--;) beta[i] = e[i+1];
      Q.DotProduct(st, z); // store the basis of deflated Krylov space temporarily in z
      for(int i = J; i--;) *Q[i] = *z[i]; // accept the new Lanczos basis
      restarted = false; // restart has been dealt with
    } // end if, we can pretend now that no restart whatsoever has ever happended
    Q.PushBackColumn(q); // confirm new Lanczos basis qj
    r.Purify(Q); // r -= Q*Q^T*r, correct roundoff error for loss of orthogonality
    beta.push_back(r.norm()); // betaj=|r|
    matrix s(jj,jj); s = 1.; // stores the eigenvetors of Tk
    TADiagonalize::TridiagQLImplicit(d,e,jj,s); // diagonalize Tk
    TADiagonalize::EigenSort(d,s,jj); // sort eigenpairs by eigenvalues in descending order
    double epsilon = 0.; // check for convergence
    for(int i = min(jj,J); i--;) epsilon = max(epsilon, fabs(beta[j]*s[j][i]));
    if(epsilon <= EPS || nm == j){
      TADiagonalize::EigenSort(d,s,jj);
      Q.DotProduct(s, z);
      if(nm == j) break; else return;
    }
    if(K-1 == j && !restarted && restartCnt < NMAXRST){
      // deflate j+1~K Ritz-pairs, and update the tridiagonal matrix //
      alpha.erase(alpha.begin()+J, alpha.end());
      beta.erase(beta.begin()+J, beta.end()-1); // beta_{k} reserved for future use
      Q.DotProduct(s, z); // z=Q*s: calculate all the Ritz vectors
      Q.EraseColumn(J, Q.size()); // remove the unwanted Ritz vectors
      for(int i = J; i--;) *Q[i] = *z[i]; // change qi for zi in Q
      r.Scale(1./beta[J], q); // compute q_{K+1} = r/|r|
      Q.PushBackColumn(q); // add q_{K+1} to Q
      H.DotProduct(q, r); // r = a*q_{K+1}
      for(int i = 0; i < J; i++){
        alpha[i] = d[i]; // the diagonal of Tk(0:J-1, 0:J-1)
        beta[i] = beta[J]*s[j][i]; // sigma_i: <zi|a|q_{K+1}>
        r.SelfAdd(-beta[i], *z[i]);
      } // end for over i
      alpha.push_back(r.DotProduct(q)); // alpha = <q_{K+1}|a|q_{K+1}>
      r.SelfAdd(-alpha[J], q); // remove q_{K+1} component, q_{K+2} constructed
      beta[J] = r.norm(); // r is clean now, beta_{K+1}=|r_{K+1}|
      j -= K-J; // one more thing, squeeze the loop for the purge
      j++; // j++: restarting itself moves a step forward in the Lanczos iteration
      restarted = true; restartCnt++;
      TAException::Info("TALanczos", "LanczosPurge: Restart #%d.", restartCnt);
    } // end if(j == K)
    if(restartCnt >= NMAXRST)
      TAException::Info("TALanczos", "LanczosPurge: Too many restarts occurred.");
  } // end for over j
  TAException::Info("TALanczos", "LanczosPurge: Too many iterations occurred.");
} // end of member function LanczosRestart

/// restart LanczosPurge with a new initial vector orthorgonal to the existing
/// invariant space, with the previous eigenpairs stored.
/// \param nPair: number of eigenpairs wanted
void TALanczos::PurgeRestart(TAHamiltonian &H, double *d, TASparseMatrix &z, int nPair){
  TASparseMatrix zz; // temporary eigenvector container for each call of LanczosPurge
  TASparseVec x; z.Clear();
  double *dd = new double[nPair*2]; // temporary eigenvalue container for each call of LanczosPurge

  while(int(z.size()) < nPair){
    InitializeVector(x, z, nPair); // assign the initial Lanczos vector
    LanczosPurge(H,dd,zz,&x, 2);
    memcpy(d+z.size(), dd, zz.size()*sizeof(double)); z.MoveBack(zz); // stash the liberated eigenpairs
  } // end while
  EigenSort(d,z,nPair);
  delete dd;
} // end of member function PurgeRestart

/// devise a vector x, such that it is orthogonal to span{z}, and has increasing
/// weight upwards (i.e., the former coordinates have more weight than the latter)
/// \param n: the dimension of the matrix to be diagonalized
void TALanczos::InitializeVector(TASparseVec &x, TASparseMatrix &z, int n){
  static TRandom3 r;
  const int nn = n < 50 ? n : 50; // no more than 50 may be a good choice

  double c = 6./((nn-1)*(2*n-1)); // rough normalization factor of x
  x.Clear();
  for(int i = 0; i < nn; i++) x.Fill(i, r.Gaus((nn-i), 0.5)*c); // Exp(tau): exp(-x/tau)
  x.Purify(z);
} // end of member function InitializeVector

/// note that eigenvalues are not ordered in Jacobi. This routine below sorts
/// the input eigenvalues into descending order, and rearranges the eigenvectors
/// correspondingly. The method is straight insertion
void TALanczos::EigenSort(double *d, TASparseMatrix &v, int n){
  double p; TASparseVec *vp; // temporary variables for data exchange
  int k; // subscript of the maximum
  for(int i = 0; i < n-1; i++){
    p = d[k=i];
    for(int j = i+1; j < n; j++) if(fabs(d[j]) >= p) p = d[k=j];
    if(k != i){
      d[k] = d[i]; d[i] = p; // exchange di and dk
      vp = v[i]; v[i] = v[k]; v[k] = vp; // exchange column i and column k
    } // end if
  } // end for over i
} // end of member function EigenSort
