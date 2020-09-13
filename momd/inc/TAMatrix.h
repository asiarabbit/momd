/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMatrix.h
  \class TAMatrix
  \brief This class just represents a matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/01
  \date Last modified: 2020/09/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMatrix_h_
#define _TAMatrix_h_

#include <cstring>
#include <initializer_list>

class TAMatrix{
public:
  TAMatrix(int nr_, int nc_);
  TAMatrix(const TAMatrix &m); // the copy constructor
  TAMatrix(TAMatrix &&m); // the move constructor
  TAMatrix &operator=(const TAMatrix &m); // the assignment constructor
  TAMatrix &operator=(TAMatrix &&m); // the assignment constructor
  TAMatrix &operator=(const std::initializer_list<double> &li);

  virtual ~TAMatrix();

  double *operator[](int i) const{ return a+i*nc; }
  bool DimensionMatch(const TAMatrix &b){ return b.nr == nr && b.nc == nc;}
  void Print();
protected:
  double *a;
  int nr, nc;
};

typedef TAMatrix matrix;

#endif
