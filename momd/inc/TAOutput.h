/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAOutput.h
  \class TAOutput
  \brief This class is responsible for printing format of the momdis results on
  screen and into file. So this is a tool class, and basically a static class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/07/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#define _TAOutput_h_
#define _TAOutput_h_

class TAOutput{
public:
  ~TAOutput(){}

  /// dedicated for printing all parts of the knockout c.s.
  static void PrintKOCS(PrintParallelMOMDIS(int fl,
    const double *sigmaStr_M, double sigmaStr,
    const double *sigmaDiff_M, double sigmaDiff,
    const double *sigmaStrTotal_M, double sigmaTotal);
  /// a file of two columns with file name being filename
  /// \param len: length of array x and y
  static void PrintToFile(int len, const double *x, const double *y,
    const string &filename);
};

#endif
