/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TASparseVec.cxx
  \class TASparseVec
  \brief Serving as a column for sparse matrix, inherited from ROOT main I/O class
   - TTree, which is dedicated to handle large amount of data. Note that one object
  stands for one column in a sparse matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/10
  \date Last modified: 2021/02/10 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include "TAMBSDTree.h"
// #include "TSQLResult.h"
// #include "TSQLRow.h"

TAMBSDTree::TAMBSDTree(const string &name, const string &title) : TATreeCol(name, title){
  // Branch("r", &fR, "r/l");
  fTree->Branch("bit", &fBit, "bit/l");
} // end of constructor}

unsigned long TAMBSDTree::GetMBSDInBit(unsigned long mbsdIndex){
  fTree->GetEntry(mbsdIndex);
  return fBit;
} // end of member function GetElement
