/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TATreeCol.h
  \class TATreeCol
  \brief storing vector elements in a TTree object for large vectors. Note that
  this is a base class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/10
  \date Last modified: 2021/02/10 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TATreeCol_h_
#define _TATreeCol_h_

#include <string>
#include "TTree.h"

using std::string;

class TFile;

class TATreeCol{
public:
  TATreeCol(const string &name, const string &title);
  TATreeCol(const TATreeCol &t); ///< the copy constructor
  virtual ~TATreeCol();
  virtual unsigned long GetEntries(){ return fTree->GetEntries(); }
  virtual void Save(); ///< save fTree
  static void CloseFile(); ///< close kCurrentFile

protected:
  static TFile *kCurrentFile; ///< rootfile to write the cached tree to
  TTree *fTree;
};

#endif
