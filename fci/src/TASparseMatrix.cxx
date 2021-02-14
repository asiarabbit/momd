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


TASparseMatrix::TASparseMatrix(){}
TASparseMatrix::~TASparseMatrix(){}

/// r = (*this)*s
void TASparseMatrix::DotProduct(const vec_t<double> &s, TASparseVec &r){
  const unsigned long m = size();
  if(s.size() != m) TAException::Error("TASparseMatrix", "Dimension Mismatch.");

  unsigned long n1 = (*std::max_element(begin(), end(),
    [](TASparseVec *a, TASparseVec *b){ return a->GetRow(0) < b->GetRow(0); }))->GetRow(0);
  unsigned long n2 = (*std::max_element(begin(), end(),
    [](TASparseVec *a, TASparseVec *b){ return a->GetMaxRow() > b->GetMaxRow(); }))->GetMaxRow();

  for(unsigned long i = n1; i < n2; i++){
    double t = 0.;
    for(unsigned long j = 0; j < m; j++) t += (*this->at(j))[i]*s[j];
    r.Fill(i, t);
  }
} // end of member function DotProduct
// r = (*this)*s
void TASparseMatrix::DotProduct(const matrix &s, TASparseMatrix &r){
  const int n = size();
  for(int i = 0; i < n; i++) DotProduct(s.cv(i), *r[i]);
}

void TASparseMatrix::PushBackColumn(TASparseVec &r){
  push_back(new TASparseVec(r));
}
/// Erase [c1, c2)
void TASparseMatrix::EraseColumn(int c1, int c2){
  std::for_each(begin()+c1, begin()+c2, [](TASparseVec *c){ delete c; });
  erase(begin()+c1, begin()+c2);
}
