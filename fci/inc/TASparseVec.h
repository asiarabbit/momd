/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TASparseVec.h
  \class TASparseVec
  \brief Serving as a column for sparse matrix, inherited from ROOT main I/O class
   - TTree, which is dedicated to handle large amount of data. Note that one object
  stands for one column in a sparse matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/10
  \date Last modified: 2021/02/12 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASparseVec_h_
#define _TASparseVec_h_

#include <string>
#include "TATreeCol.h"
#include "TABit.h"
#include "TTree.h"

class TASparseMatrix;

using std::string;

class TASparseVec : public TATreeCol{
public:
  TASparseVec(const string &name = "v", const string &title = "");
  TASparseVec(const TASparseVec &v);
  virtual ~TASparseVec(){}

  double GetValue(unsigned long r);
  double operator[](unsigned long r){ return GetValue(r); }
  unsigned long GetRow(unsigned long entryIndex); ///< \retval fR of entryIndex
  void Clear();
  void Fill(unsigned long r, double v);
  void SelfAdd(double k, TASparseVec &v); ///<\retval r += k*v
  void Scale(double v, TASparseVec &r); ///< r = v*this
  double DotProduct(TASparseVec &v); ///< \retval (*this)*v
  void Purify(TASparseMatrix &Q); ///< remove Q=(q1,q2..qn) bits from r: r-= \sum_q{|q><q|r>}
  unsigned long GetMinRow(); ///< \retval row of the first (non-zero) entry
  unsigned long GetMaxRow(); ///< \retval row of the last (non-zero) entry
  double norm(); ///< \retval sqrt((*this)*(*this))
  void operator=(const TASparseVec &v); ///< copy v to this

protected:
  unsigned long fR; ///< row number
  double fV; ///< the matrix element
  char fSelection[512]; ///< ROOT selection, e.g. "r==3"
};

#endif
