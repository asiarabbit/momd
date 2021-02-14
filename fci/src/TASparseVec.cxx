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
#include "TASparseVec.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TASparseMatrix.h"

TASparseVec::TASparseVec(const string &name, const string &title) : TATreeCol(name, title){
  fTree->Branch("r", &fR, "r/l");
  fTree->Branch("v", &fV, "v/D");
} // end of constructor

TASparseVec::TASparseVec(const TASparseVec &t) : TATreeCol(t){}

void TASparseVec::Fill(unsigned long r, double v){
  if(0. == v) return;
  fR = r; fV = v;
  fTree->Fill();
}

double TASparseVec::GetValue(unsigned long r){
  sprintf(fSelection, "r==%lu", r);
  return std::stod(fTree->Query("v", fSelection)->Next()->GetField(0)); // 7 digits in decimal part
} // end of member function GetElement

void TASparseVec::Clear(){
  TTree *t = fTree->CloneTree(0); ///< 0: clone all structure (branches) but no entries
  delete fTree; fTree = t;
}

/// \retval fR of entryIndex
unsigned long TASparseVec::GetRow(unsigned long entryIndex){
  fTree->GetEntry(entryIndex);
  return fR;
}
///\retval r += k*v
void TASparseVec::SelfAdd(double k, TASparseVec &v){
  // obtain the entry index range with fV != 0. or v.fV != 0. //
  const unsigned long n1 = std::min(GetRow(0), v.GetRow(0));
  const unsigned long n2 = std::max(GetMaxRow(), v.GetMaxRow());

  TTree *t = fTree->CloneTree(0); ///< 0: clone all structure (branches) but no entries
  for(unsigned long i = n1; i <= n2; i++){
    fR = i; fV = this->GetValue(i) + k*v[i];
    t->Fill();
  }
  delete fTree; fTree = t;
} // end member function SelfAdd
/// r = k*this
void TASparseVec::Scale(double k, TASparseVec &r){
  r.Clear(); // remove all the entries in r
  unsigned long n = fTree->GetEntries();
  for(unsigned long i = 0; i < n; i++){
    fTree->GetEntry(i);
    r.Fill(fR, k*fV);
  }
} // end member function Scale

/// \retval (*this)*v
double TASparseVec::DotProduct(TASparseVec &v){
  const unsigned long n1 = std::min(GetRow(0), v.GetRow(0));
  const unsigned long n2 = std::max(GetMaxRow(), v.GetMaxRow());
  double r = 0., e1, e2;
  for(unsigned long i = n1; i <= n2; i++)
    if((e1 = this->GetValue(i)) && (e2 = v[i])) r += e1 * e2;
  return r;
}
/// \retval (*this)*(*this)
double TASparseVec::norm(){
  const unsigned long n1 = GetRow(0);
  const unsigned long n2 = GetMaxRow();
  double r = 0., e1;
  for(unsigned long i = n1; i <= n2; i++)
    if((e1 = this->GetValue(i))) r += e1 * e1;
  return sqrt(r);
}

/// remove Q=(q1,q2..qn) bits from r: r-= \sum_q{|q><q|r>}
void TASparseVec::Purify(TASparseMatrix &Q){
  for(TASparseVec *q : Q) SelfAdd(-DotProduct(*q), *q);
}

/// \retval row of the last (non-zero) entry
unsigned long TASparseVec::GetMaxRow(){ return GetRow(fTree->GetEntries()-1); }

/// copy v to this
void TASparseVec::operator=(const TASparseVec &v){
  TTree *t = (TTree *)v.fTree->Clone("assign");
  delete fTree; fTree = t;
}
