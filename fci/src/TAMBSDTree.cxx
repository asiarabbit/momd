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
  CreateBranch();
} // end of constructor}

unsigned long long TAMBSDTree::GetMBSDInBit(unsigned long long mbsdIndex){
  GetEntry(mbsdIndex);
  return fBit;
} // end of member function GetElement

void TAMBSDTree::Fill(unsigned long long bit){
  SetBranchAddress();
  fBit = bit; fTree->Fill();
}

void TAMBSDTree::CreateBranch(){
  fTree->Branch("bit", &fBit, "bit/l");
}
void TAMBSDTree::SetBranchAddress(){
  if((void*)(&fBit) == (void*)fTree->GetBranch("bit")->GetAddress()) return; // already set

  fTree->SetBranchAddress("bit", &fBit);
}
