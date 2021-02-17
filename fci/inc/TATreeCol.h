/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TATreeCol.h
  \class TATreeCol
  \brief storing vector elements in a TTree object for large vectors. Supposed
  to be a base class.
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
  /// NOTE that user has to reset branch addresses to their own
  // member data manually in the derived classes
  TATreeCol &operator=(const TATreeCol &t);
  TATreeCol(const TATreeCol &t);
  virtual ~TATreeCol();
  virtual unsigned long long GetEntries() const{ return fTree->GetEntries(); }
  virtual void Save(); ///< save fTree
  static void CloseFile(); ///< close kCurrentFile
  virtual void SetBranchAddress() = 0; ///< Set Branch address to fTree
  virtual void CreateBranch() = 0; ///< create branches
  virtual void GetEntry(unsigned long long i); ///< TTree::GetEntry with branch address checking

protected:
  static TFile *kCurrentFile; ///< rootfile to write the cached tree to
  TTree *fTree;
  static int kObjCount; ///< to give each TATreeCol obj a globally unique id
};

#endif
