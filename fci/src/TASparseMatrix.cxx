/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TASparseMatrix.cxx
  \class TASparseMatrix
  \brief dedicated for Storage of sparse matrix. The class is a vector of matrix
  columns, with each column stored in a ROOT TTree object - TATreeCol
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/12
  \date Last modified: 2021/02/12 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <algorithm>
#include "TASparseMatrix.h"
#include "TASparseVec.h"
#include "TAException.h"

using std::for_each;

TASparseMatrix::TASparseMatrix(){}
TASparseMatrix::~TASparseMatrix(){}

/// r = (*this)*s
void TASparseMatrix::DotProduct(const vec_t<double> &s, TASparseVec &r){
  const unsigned long long m = size();
  if(s.size() != m) TAException::Error("TASparseMatrix", "Dimension Mismatch.");

  unsigned long long n1 = (*std::min_element(begin(), end(),
    [](TASparseVec *a, TASparseVec *b){ return a->GetRow(0) < b->GetRow(0); }))->GetRow(0);
  unsigned long long n2 = (*std::max_element(begin(), end(),
    [](TASparseVec *a, TASparseVec *b){ return a->GetMaxRow() < b->GetMaxRow(); }))->GetMaxRow();

  double e1, e2;
  for(unsigned long long i = n1; i <= n2; i++){
    double t = 0.;
    for(unsigned long long j = 0; j < m; j++) if((e1 = s[j]) && (e2 = (*this->at(j))[i])) t += e1*e2;
    if(t) r.Fill(i, t);
  }
} // end of member function DotProduct
// r = (*this)*s
void TASparseMatrix::DotProduct(const matrix &s, TASparseMatrix &r){
  if(!r.size()){
    r.Clear(); // empty r
    // TAException::Error("TASparseMatrix", "DotProduct: r is supposed to be empty");
  }
  const int n = size();
  for(int i = 0; i < n; i++){
    r.push_back(new TASparseVec("r"));
    DotProduct(s.cv(i), *r[i]);
  }
} // end member function DotProduct

void TASparseMatrix::PushBackColumn(const TASparseVec &r){
  push_back(new TASparseVec(r));
}
/// Erase [c1, c2)
void TASparseMatrix::EraseColumn(int c1, int c2){
  if(c2 >= int(size()) || c1 < 0)
    TAException::Error("TASparseMatrix", "EraseColumn: Designated range out of border.");
  for_each(begin()+c1, begin()+c2, [](TASparseVec *c){ delete c; });
  erase(begin()+c1, begin()+c2);
}

void TASparseMatrix::Print() const{
  for_each(begin(), end(), [](TASparseVec *c){ c->Print(); });
}

void TASparseMatrix::Clear(){
  if(!size()) return;
  EraseColumn(0, size());
}

/// only transfer ownership. r is emptied.
void TASparseMatrix::MoveBack(TASparseMatrix &r){
  for_each(r.begin(), r.end(), [this](TASparseVec *c){ push_back(c); });
  r.clear();
}

void TASparseMatrix::Save(){
  for_each(begin(), end(), [](TASparseVec *c){ c->Save(); });
}
