/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMOMDIS_M.h
  \class TAMOMDIS_M
  \brief Calculate residue momentum distribution with angular momentum component
  m specified. This is a class to assist class TAMOMDIS.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/09/04 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMOMDIS_M_h_
#define _TAMOMDIS_M_h_

class TAMOMDIS;

class TAMOMDIS_M{
public:
  TAMOMDIS_M(TAMOMDIS *mom);
  virtual ~TAMOMDIS_M();

  /// returning c.s. for m, mom is assigned to momStr
  double ParallelStr(int m, double *momStr);
  /// the same as ParallelStr, but for diffraction dissociation
  double ParallelDiff(int m, double *momDiff);

  double *GetMomArr() const{ return fMomArr; }

  static const int kNmom = 200; ///< number of points in momentum distribution

private:
  TAMOMDIS *fMOM;
  double *fMomArr; ///< momentum array for the momentum distribution
};

#endif
