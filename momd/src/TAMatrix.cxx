/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMatrix.cxx
  \class TAMatrix
  \brief This class just represents a matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/01
  \date Last modified: 2020/09/08 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <cstdio>
#include "TAMatrix.h"
#include "TAException.h"

TAMatrix::TAMatrix(int nr_, int nc_) : nr(nr_), nc(nc_){
  if(!nr || !nc) a = nullptr;
  else a = new double[nr*nc]{};
} // end of the constructor
TAMatrix::~TAMatrix(){
  if(a) { delete [] a; a = nullptr; }
} // end of the destructor

TAMatrix &TAMatrix::operator=(const TAMatrix &m){
  if(this == &m) return *this;
  if(m.nr != this->nr || m.nc != this->nc)
    TAException::Error("TAMatrix", "assignment_ctor: dimension mismatch.");

  memcpy(a, m.a, nr*nc*sizeof(double));
  return *this;
} // end of the assignment constructor
TAMatrix &TAMatrix::operator=(const std::initializer_list<double> &li){
  const int n = nr * nc;
  if(int(li.size()) != n)
    TAException::Error("TAMatrix", "operator={}: initializer_list size\
doesn't match the number of the matrix element.");

  int cnt = 0;
  for(auto t : li){
    a[cnt++] = t;
    if(cnt >= n) break;
  } // end for
  return *this;
} // end assignment with {}

TAMatrix::TAMatrix(const TAMatrix &m){
  *this = m;
} // end of the copy constructor

TAMatrix::TAMatrix(TAMatrix &&m){
  *this = m;
} // end of the move constructor
TAMatrix &TAMatrix::operator=(TAMatrix &&m){
  if(m.nr != this->nr || m.nc != this->nc)
    TAException::Error("TAMatrix", "move_assignment: dimension mismatch.");

  if(a) delete [] a;
  a = m.a; m.a = nullptr;
  return *this;
} // end of the move constructor

void TAMatrix::Print(){
  for(int i = 0; i < nr; i++){
    for(int j = 0; j < nc; j++){
      printf(" %10.6f", (*this)[i][j]);
    } // end for over columns
    printf("\n");
  } // end loop over rows
  printf("\n");
} // end of member function Print
