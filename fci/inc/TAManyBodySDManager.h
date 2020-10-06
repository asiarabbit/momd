/**
  \file TAManyBodySDManager.h
  \class TAManyBodySDManager
  \brief A class to generate many-body basis and manage TAManagerSDList objects.
  \author SUN Yazhou
  \date Created: 2020/02/01
  \date Last modified: 2020/02/10 by SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAManyBodySDManager_h_
#define _TAManyBodySDManager_h_

#include <vector>
#include <list>

using std::vector;
using std::list;

class TAManyBodySD;
class TAManyBodySDList;

class TAManyBodySDManager{
public:
  virtual ~TAManyBodySDManager();
  static TAManyBodySDManager *Instance();
  void GenerateManyBodySD();
  void MSchemeGo(); ///< generate the M-scheme many-body state basis
  TAManyBodySDList *GetMBSDListM();

protected:
  TAManyBodySDManager();

  static TAManyBodySDManager *kInstance;
  vector<TAManyBodySD *> fManyBodySDVec; ///< the total MBSDs
  TAManyBodySDList *fManyBodySDListM; ///< M-scheme many-body basis
};

#endif
