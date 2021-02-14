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

TATreeCol::TATreeCol(const string &name, const string &title){
  if(!kCurrentFile) kCurrentFile = new TFile("fci.root", "RECREATE");
  fTree = new TTree(name.c_str(), title.c_str());
} // end of constructor

/// the copy constructor
TATreeCol::TATreeCol(const TATreeCol &t){
  if(!kCurrentFile) kCurrentFile = new TFile("fci.root", "UPDATE");
  fTree = (TTree *)t.fTree->Clone("clone");
}

TATreeCol::~TATreeCol(){
  if(fTree){
    delete fTree; fTree = nullptr;
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
