/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAIntegral.cxx
  \class TAIntegral
  \brief This class implements general numerical integral. Input and results are
  all stored in arrays. Note that this is a mathematical tool class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/09/03 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include "TAIntegral.h"

// integral of f(x) over domain [x[0],x[n-1]]. n is the length of x and f
// composed with formula (8.3.6) p.212 Computing Method ver.3 by Guicheng Li
template<typename T, typename FUNCTOR>
T TAIntegral<T, FUNCTOR>::Simpson(int n, const double *x, const T *f){
  // number of individual Simpson intervals
  int k = (n - 1) / 2; // f[n-1] is dropped in the case where n is even

  T sum = f[0] + f[2*k];
  for(int i = 0; i < k; i++) sum += 4.*f[2*i+1];
  for(int i = 1; i < k; i++) sum += 2.*f[2*i];

  const double h = (x[n-1] - x[0]) / (n - 1);
  if(h < 0){
    printf("Array x should be in ascending order!!!\n");
    printf("Exiting...\n");
    exit(1);
  }
  sum *=  h / 3.;
  // Simpson's rule is only for odd n
  // otherwise integral of f[n-1]*h should be explicitly included
  if(n % 2 == 0) sum += f[n-1] * h;

  return sum;
} // end of member function Simpson
