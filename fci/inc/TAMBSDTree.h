/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TAMBSDTree.h
  \class TAMBSDTree
  \brief Serving as a vector to store many-body Slater determinants.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/10
  \date Last modified: 2021/02/10 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMBSDTree_h_
#define _TAMBSDTree_h_

#include <string>
#include "TATreeCol.h"
#include "TABit.h"
#include "TTree.h"

using std::string;

class TAMBSDTree : public TATreeCol{
public:
  TAMBSDTree(const string &name = "mbsdList", const string &title = "MBSDList");
  virtual ~TAMBSDTree(){}

  unsigned long GetMBSDInBit(unsigned long mbsdIndex);
  unsigned long operator[](unsigned long mbsdIndex){ return GetMBSDInBit(mbsdIndex); }
  void Fill(unsigned long bit){ fBit = bit; fTree->Fill(); }

protected:
  // unsigned long fR; ///< row number
  unsigned long fBit; ///< the many-body Slater determinant in bit representation
  // char fSelection[512]; ///< ROOT selection, e.g. "r==3"
};

#endif
