/**
  \file TAManyBodySDManager.h
  \class TAManyBodySDManager
  \brief A class to generate many-body basis and manage TAManyBodySDList objects.
  \author SUN Yazhou
  \date Created: 2020/02/01
  \date Last modified: 2020/02/10 by SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAManyBodySDManager_h_
#define _TAManyBodySDManager_h_

#include <vector>
#include <string>

using std::vector;
using std::string;

class TAManyBodySD;
class TAManyBodySDList;

class TAManyBodySDManager{
public:
  virtual ~TAManyBodySDManager();
  static TAManyBodySDManager *Instance();
  /// assign nparticle, nspstate, 2*M and the spstate space
  void LoadConfigFile(const string &file);
  void GenerateManyBodySD();
  void MSchemeGo(); ///< generate the M-scheme many-body state basis
  TAManyBodySDList *GetMBSDListM();
  int GetNParticle() const{ return fNParticle; }
  int GetNSPState() const{ return fNSPState; }
  int Get2M() const{ return f2M; }

protected:
  TAManyBodySDManager();

  static TAManyBodySDManager *kInstance;
  TAManyBodySDList *fManyBodySDListM; ///< M-scheme many-body basis
  int fNParticle, fNSPState, f2M;
  string fSPStatefile;
};

#endif
