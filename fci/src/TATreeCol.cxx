/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TATreeCol.cxx
  \class TATreeCol
  \brief storing vector elements in a TTree object for large vectors
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2021/02/10
  \date Last modified: 2021/02/10 by SUN Yazhou
  \copyright 2020-2021 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include "TFile.h"
#include "TATreeCol.h"
#include "TAException.h"

TFile *TATreeCol::kCurrentFile = nullptr;
int TATreeCol::kObjCount = 0;

TATreeCol::TATreeCol(const string &name, const string &title){
  if(!kCurrentFile) kCurrentFile = new TFile("fci.root", "RECREATE");
  char name1[512];
  sprintf(name1, "%s_%d", name.c_str(), kObjCount++);
  fTree = new TTree(name1, title.c_str());
} // end of constructor

/// NOTE that user has to set branch addresses to their own member data manually
TATreeCol::TATreeCol(const TATreeCol &t){
  fTree = nullptr;
  *this = t;
}
TATreeCol &TATreeCol::operator=(const TATreeCol &t){
  TTree *tt = t.fTree->CloneTree(t.GetEntries());
  if(fTree) delete fTree;
  fTree = tt;
  return *this;
}

TATreeCol::~TATreeCol(){
  if(fTree){
    delete fTree;
    fTree = nullptr;
  }
} // end of destructor

void TATreeCol::Save(){
  fTree->Write("", TObject::kOverwrite);
}

void TATreeCol::CloseFile(){
  if(!kCurrentFile) TAException::Error("TATreeCol", "CloseFile: Called once already.");
  kCurrentFile->Close();
  delete kCurrentFile; kCurrentFile = nullptr;
}

/// TTree::GetEntry with branch address checking
void TATreeCol::GetEntry(unsigned long long i){
  SetBranchAddress();
  if(i > GetEntries()) TAException::Error("TATreeCol", "GetEntry: required entry index out of range");
  fTree->GetEntry(i);
}
