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

using std::vector;
using std::min;
using std::max;

/// Lanczos algorithm is devised for a few largest (in modulus) eigenvalues and
/// It projects A to a space with much lower dimension (Krylov subspace) as a tridiagonal
/// matrix, whose eigenvalues and eigenvectors are good approximation of A's
/// Given input matrix a[0..n-1][0..n-1], the routine returns first fiew (usually <=10)
/// eigenvalues in d and the corresponding eigenvectors in z
/// d should be at least of length 50, for it is drafted in the program
void TALanczos::Lanczos(TAHamiltonian &H, double *d, TASparseMatrix &z){
  static const double EPS = 1.e-6; // eigenvalue accuracy standard
  static const int NMAX = 50; // maximum number of iterations
  const int nm = NMAX, np = nm + 1;
  double e[np]; vector<double> alpha, beta;
  // prepare for the start of the iterations //
  TASparseVec v("v"), q("q"), r("r"); // q_{j-1}, qj, q_{j+1}, 3 basis vectors
  q.Fill(0, 1.); // initialize q to (1,0,0,..)
  TASparseMatrix Q; Q.PushBackColumn(q); // to store the Lanczos vectors
  H.DotProduct(q, r); // r = a*q
  alpha.push_back(r.DotProduct(q)); // alpha = <qj|a|qj>
  r.SelfAdd(-alpha[0], q); // r -= alphaj*qj
  beta.push_back(r.norm()); // the first off-diagonal of Tk
  // the trial vector is an eigenvector
  if(fabs(beta[0])+1. == 1.){ d[0] = alpha[0]; return; }
  // OK, here we commence the iteration //
  for(int j = 1; j <= nm; j++){
    v = q; // q_{j-1}
    r.Scale(1./beta[j-1],q); // normalize qj in q
    Q.PushBackColumn(q);
    H.DotProduct(q, r); // r = a*qj
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
      if(nm == j) break; else return;
    }
  } // end for over j
  TAException::Info("TADiagonalize", "Lanczos: Too many iterations occurred.");
} // end of member function Lanczos

/// Lanczos method featuring Ritz-pairs purging. After getting a set of Ritz vectors
/// from a plain Lanczos, we keep a smaller set and remove (purge) the rest. Then
/// develop Krylov space from the new Lanczos basis. So unwanted directions are removed
/// from the initial vector for each iteration. NOTE that the largest eigenvalues are kept.
/// Ref.: http://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf, p212
/// Ref.: K. Wu and H. D. Simon,SIAM J. Matrix Anal. Appl., 22 (2000), pp. 602–616.
void TALanczos::LanczosPurge(TAHamiltonian &H, double *d, TASparseMatrix &z){
  static const double EPS = 1.e-8; // eigenvalue accuracy standard
  static const int NMAXRST = 100; // maximum restart times
  // number of wanted eigenpairs and dimension for the plain Lanczos
  static const int J = 4, K = J + 3; // J, K: number of the wanted eigenpairs and candidates
  const unsigned long n = H.GetNBasis();
  if(n <= 7) return Lanczos(H,d,z); // we deal with a little bit larger matrix
  const int NMAX = 100; // maximum number of iterations
  const int nm = min(n, (unsigned long)NMAX);
  const int np = nm+1; // maximum dimension of the final Tk
  double e[np]{}; vector<double> alpha, beta;
  int restartCnt = 0; // restart times
  bool restarted = false; // if has just restarted
  // prepare for the start of the iterations //
  matrix st(J+1,J+1); st = 1.; // s: the eigenvector matrix of Tk, st set to E
  TASparseVec v("v"), q("q"), r("r"); // q_{j-1}, qj, q_{j+1}, 3 basis vectors
  q.Fill(0, 1.); // initialize q to (1,0,0,..)
  TASparseMatrix Q; Q.PushBackColumn(q); // to store the Lanczos vectors
  H.DotProduct(q, r); // r = a*q
  alpha.push_back(r.DotProduct(q)); // alpha = <qj|a|qj>
  r.SelfAdd(-alpha[0], q); // r = A*q-alphaj*qj
  beta.push_back(r.norm()); // the first off-diagonal of Tk
  // FIXME: if the trial vector is an eigenvector, return
  if(fabs(beta[0])+alpha[0] == alpha[0]){ d[0] = alpha[0]; return; }
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
      // cout << "epsilon: \033[34;1m" << epsilon << "\033[0m\n"; // DEBUG
      // printf("restartCnt: %d, j = %d\n", restartCnt, j); // DEBUG
      // for(int i = 0; i < jj; i++) printf("d[%d]: %f\n", i, d[i]); // DEBUG
      // for(int i = 0; i < jj; i++) printf("beta[%d]*s[%d][%d]: %f\n", j,j,i,beta[j]*s[j][i]); // DEBUG
      // getchar(); // DEBUG
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
    } // end if(j == K)
  } // end for over j
  TAException::Info("TADiagonalize", "LanczosPurge: Too many iterations occurred.");
} // end of member function LanczosRestart
