/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TASparseVec.cxx
  \class TASparseVec
  \brief Serving as a column for sparse matrix, inherited from ROOT main I/O class
   - TTree, which is dedicated to handle large amount of data. Note that one object
  stands for one column in a sparse matrix, so each entry is a 3-member tuple. Zero
  elements are just skipped.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/10
  \date Last modified: 2021/02/10 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <algorithm>
#include <string>
#include <cmath>
#include <iomanip>
#include "TASparseVec.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TASparseMatrix.h"

using std::cout;
using std::endl;
using std::setw;
using std::ios_base;

TASparseVec::TASparseVec(const string &name, const string &title) : TATreeCol(name, title){
  CreateBranch();
} // end of constructor

TASparseVec::TASparseVec(const TASparseVec &t) : TATreeCol(t){
  SetBranchAddress();
}
/// copy v to this
TASparseVec &TASparseVec::operator=(const TASparseVec &v){
  TATreeCol::operator=(v);
  SetBranchAddress();
  return *this;
} // end of assignment function


void TASparseVec::CreateBranch(){
  fTree->Branch("index", &fIndex, "index/l");
  fTree->Branch("r", &fR, "r/l");
  fTree->Branch("v", &fV, "v/D");
}
void TASparseVec::SetBranchAddress(){
  if((void*)(&fV) == (void*)fTree->GetBranch("v")->GetAddress()) return; // already set

  fTree->SetBranchAddress("index", &fIndex);
  fTree->SetBranchAddress("r", &fR);
  fTree->SetBranchAddress("v", &fV);
}

/// \retval row of the last (non-zero) entry
unsigned long long TASparseVec::GetMaxRow(){
  return GetRow(fTree->GetEntries()-1);
}

void TASparseVec::Fill(unsigned long long r, double v){
  if(0. == v) return;

  fIndex = fTree->GetEntries();
  fR = r; fV = v;
  SetBranchAddress();
  fTree->Fill();
}

double TASparseVec::GetValue(unsigned long long r){
  char selection[512];
  sprintf(selection, "r==%llu", r);
  TSQLResult *res = fTree->Query("index", selection);
  if(!res->GetRowCount()) return 0.; // found none
  GetEntry(std::stoul(res->Next()->GetField(0))); delete res;
  return fV;
} // end of member function GetElement

void TASparseVec::Clear(){
  TTree *t = fTree;
  fTree = new TTree(t->GetName(), t->GetTitle());
  t->Delete("all");

  CreateBranch();
}

/// \retval fR of entryIndex
unsigned long long TASparseVec::GetRow(unsigned long long entryIndex){
  if(entryIndex >= (unsigned long long)fTree->GetEntries()){
    TAException::Warn("TASparseVec", "GetRow: Required entry index out of border.");
    return 0;
  }

  GetEntry(entryIndex);
  return fR;
}
///\retval r += k*v
void TASparseVec::SelfAdd(double k, TASparseVec &v){
  // obtain the entry index range with fV != 0. or v.fV != 0. //
  if(!v.GetEntries() || !k) return;
  if(!GetEntries()) v.Scale(k, *this); // r = 0, and r += k*v => r = k*v
  const unsigned long long n1 = std::min(GetRow(0), v.GetRow(0));
  const unsigned long long n2 = std::max(GetMaxRow(), v.GetMaxRow());

  TTree *t = fTree->CloneTree(0);
  for(unsigned long long i = n1; i <= n2; i++){
    fV = GetValue(i)+k*v[i];
    if(fV){ fIndex = t->GetEntries(); fR = i; t->Fill(); }
  }
  delete fTree; fTree = t;
  SetBranchAddress();
  // fTree->Scan(); Print(); // DEBUG
} // end member function SelfAdd
/// r = k*this
void TASparseVec::Scale(double k, TASparseVec &r){
  r.Clear(); // remove all the entries in r
  if(!GetEntries()) return; // this is a zero vector, so just leave r empty is enough
  const unsigned long long n = fTree->GetEntries();
  for(unsigned long long i = 0; i < n; i++){
    GetEntry(i);
    r.Fill(fR, k*fV);
  }
} // end member function Scale

/// \retval (*this)*v
double TASparseVec::DotProduct(TASparseVec &v){
  if(!GetEntries()) return 0.;
  const unsigned long long n1 = std::min(GetRow(0), v.GetRow(0));
  const unsigned long long n2 = std::max(GetMaxRow(), v.GetMaxRow());
  double r = 0., e1, e2;
  for(unsigned long long i = n1; i <= n2; i++)
    if((e1 = GetValue(i)) && (e2 = v[i])) r += e1 * e2;
  return r;
}
/// \retval (*this)*(*this)
double TASparseVec::norm(){
  if(!GetEntries()) return 0.;
  const unsigned long long n1 = GetRow(0);
  const unsigned long long n2 = GetMaxRow();
  double r = 0., e1;
  for(unsigned long long i = n1; i <= n2; i++) if((e1 = GetValue(i))) r += e1 * e1;
  return sqrt(r);
}

/// remove Q=(q1,q2..qn) bits from r: r-= \sum_q{|q><q|r>}
void TASparseVec::Purify(TASparseMatrix &Q){
  if(!GetEntries()) return; // zero vector is orthogonal to any vector, need no reorthogonalizations
  for(TASparseVec *q : Q) SelfAdd(-DotProduct(*q), *q);
}

void TASparseVec::Print(){
  // cout << "Branch status - index: " << fTree->GetBranchStatus("index") << endl;
  // cout << "Branch status - r: " << fTree->GetBranchStatus("r") << endl;
  // cout << "Branch status - v: " << fTree->GetBranchStatus("v") << endl;
  const int n = 20 < GetEntries() ? 20 : GetEntries();
  cout << "TASparseVec Print: totally " << n << " elements." << endl;
  ios_base::fmtflags initial = cout.setf(ios_base::fixed, ios_base::floatfield);
  cout.unsetf(ios_base::floatfield);
  cout.precision(6);
  cout << std::right;
  for(int i = 0; i < n; i++){
    GetEntry(i);
    cout << setw(10) << fR;
  }
  cout << "\033[32;1m" << endl;
  for(int i = 0; i < n; i++){
    GetEntry(i);
    cout << setw(10) << fV;
  }
  cout << "\033[0m" << endl;
  cout.setf(initial);
  // fTree->Print();
  // fTree->Scan();
} // end of member function Print

void TASparseVec::normalize(){
  if(!GetEntries()) return;
  const double p = norm();
  if(1 == p) return;

  TASparseVec tmp;
  Scale(1./p, tmp);
  delete fTree;
  fTree = tmp.fTree;
  tmp.fTree = nullptr;
  SetBranchAddress();
}
