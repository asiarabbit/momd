/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAOutput.cxx
  \class TAOutput
  \brief This class is responsible for printing format of the momdis results on
  screen and into file. So this is a tool class, and basically a static class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/09/06 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include "TAOutput.h"

using std::cout;
using std::endl;
using std::setw;
using std::ofstream;

void TAOutput::PrintKOCS(int fl,
  const double *sigmaStr_M, double sigmaStr,
  const double *sigmaDiff_M, double sigmaDiff,
  const double *sigmaTotal_M, double sigmaTotal){
  // stripping cross secton
  cout << "Stripping cross section (mb):" << endl;
  for(int m = 0; m <= fl; m++)
    cout << setw(8) << "M: " << setw(5) << m << sigmaStr_M[m] << endl;
  cout << setw(8) << "Total: " << sigmaStr << endl;
  // diffraction cross section
  cout << "Diffraction cross section (mb):" << endl;
  for(int m = 0; m <= fl; m++)
    cout << "M: " << setw(5) << m << sigmaDiff_M[m] << endl;
  cout << setw(8) << "Total: " << sigmaDiff << endl;
  // total cross section
  cout << "Total (stripping+diffraction) cross section (mb):" << endl;
  for(int m = 0; m <= fl; m++)
    cout << "M: " << setw(5) << m << sigmaTotal_M[m] << endl;
  cout << setw(8) << "Total: " << sigmaTotal << endl;
} // end member function PrintKOCS

/// a file of two columns with file name being filename
void TAOutput::PrintToFile(int len, const double *x, const double *y,
    const string &filename){
  ofstream fout(filename.c_str());
  for(int i = 0; i < len; i++) fout << setw(10) << x[i] << setw(10) << y[i] << endl;
  fout.close();
} // end member function PrintToFile
