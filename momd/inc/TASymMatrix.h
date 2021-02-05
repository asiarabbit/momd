/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TASymMatrix.h
  \class TASymMatrix
  \brief A template class for symmetric matrix. The symmetric counterpart of an
  element is not stored. This class is only supposed to facilitate the retrieval
  of the matrix elements. Note that only the LOWER triangle is stored.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/11/25
  \date Last modified: 2020/11/25 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASymMatrix_h_
#define _TASymMatrix_h_

template <typename T>
class TASymMatrix{
public:
  TASymMatrix(int n) : fLength(n*(n+1)/2){ fData = new T[fLength]; }
  virtual ~TASymMatrix(){ delete [] fData; }

  T* operator[](int i){ return fData+i*(i+1)/2; }
  void Initialize(const T &v){ for(int i = fLength; i--;) fData[i] = v; }

private:
  const int fLength;
  T *fData;
};

#endif
