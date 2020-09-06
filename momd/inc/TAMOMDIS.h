/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMOMDIS.h
  \class TAMOMDIS
  \brief This is a global class, to take control of the whole flow of the
  program. It is also responsible for generating various momentum distributions.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/08
  \date Last modified: 2020/09/04 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#include <string>
#include <vector>

using std::string;
using std::vector;

#ifndef _TAMOMDIS_h_
#define _TAMOMDIS_h_

class TABound;
class TASMatrix;
class TAMOMDIS_M;

class TAMOMDIS{
public:
  TAMOMDIS(const string &configFile); ///< the constructor
  virtual ~TAMOMDIS();

  void LoadConfigFile(); // read configFile to fVecConfig
  void Configure(); ///< solve Rl, generate and fit S-matrix
  void Parallel(); ///< calculate dsigma/dkz
  void Go(); ///< here we go

  friend class TAMOMDIS_M;

protected:
  /// specify the config file folder
  string fConfigFile;
  vector<double> fVecConfig; // to store the content of the configFile
  /// the component to calculate momtum distribution of the core
  TABound *fBound;
  TASMatrix *fSMatrix;
  TAMOMDIS_M *fMOM; // calculate m-specific c.s.
};

#endif
