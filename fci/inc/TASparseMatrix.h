/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TASparseMatrix.h
  \class TASparseMatrix
  \brief dedicated for Storage of sparse matrix. The class is a vector of matrix
  columns, with each column stored in a ROOT TTree object - TATreeCol
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/10
  \date Last modified: 2021/02/10 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASparseMatrix_h_
#define _TASparseMatrix_h_

#include "TTree.h"
#include "TAMatrix.h"
#include <vector>

class TASparseVec;
using std::vector;

class TASparseMatrix : public vector<TASparseVec *>{
public:
  TASparseMatrix();
  virtual ~TASparseMatrix();

  void DotProduct(const vec_t<double> &s, TASparseVec &r); // r = (*this)*s
  void DotProduct(const matrix &s, TASparseMatrix &r); // r = (*this)*s
  void PushBackColumn(const TASparseVec &r);
  void MoveBack(TASparseMatrix &r); ///< only transfer ownership. r is emptied.
  void EraseColumn(int c1, int c2); ///< Erase [c1, c2)
  void Print() const;
  void Clear();
  void Save(); ///< save all the trees
};

#endif
