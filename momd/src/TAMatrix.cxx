/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMatrix.cxx
  \class TAMatrix
  \brief This class just represents a matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/01
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TAMatrix.h"
#include "TAException.h"

TAMatrix::TAMatrix(int nr_, int nc_) : nr(nr_), nc(nc_){
  a = new double[nr*nc]{};
} // end of the constructor
TAMatrix &TAMatrix::operator=(const TAMatrix &m){
  if(this == &m) return *this;
  if(m.nr != this->nr || m.nc != this->nc)
    TAException::Error("TAMatrix", "assignment_ctor: dimension mismatch.");

  memcpy(this, m.a, nr*nc*sizeof(double));
  return *this;
} // end of the assignment constructor

TAMatrix::TAMatrix(const TAMatrix &m){
  *this = m;
} // end of the copy constructor

TAMatrix::TAMatrix(TAMatrix &&m){
  if(m.nr != this->nr || m.nc != this->nc)
    TAException::Error("TAMatrix", "move_ctor: dimension mismatch.");

  delete [] this->a;
  this->a = m.a;
} // end of the move constructor
