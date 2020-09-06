/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMatrix.h
  \class TAMatrix
  \brief This class just represents a matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/08/01
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMatrix_h_
#define _TAMatrix_h_

#include <cstring>

class TAMatrix{
public:
  TAMatrix(int nr_, int nc_);
  TAMatrix(const TAMatrix &m); // the copy constructor
  TAMatrix(TAMatrix &&m); // the move constructor
  TAMatrix &operator=(const TAMatrix &m); // the assignment constructor

  virtual ~TAMatrix(){
    delete [] a; a = nullptr;
  }

  double *operator[](int i) const{ return a+i*nc; }
protected:
  double *a;
  int nr, nc;
};

typedef TAMatrix matrix;

#endif
