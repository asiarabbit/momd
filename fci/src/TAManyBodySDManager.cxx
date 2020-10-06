/**
  \file TAManyBodySDManager.C
  \class TAManyBodySDManager
  \brief A class to generate many-body basis and manage TAManagerSDList objects.
  \author SUN Yazhou
  \date Created: 2020/02/01
  \date Last modified: 2020/02/02 by SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <string>
#include <algorithm>
#include "TAManyBodySDManager.h"
#include "TAManyBodySDList.h"
#include "TAManyBodySD.h"
#include "TASingleParticleStateManager.h"
#include "TAException.h"
#include "TAMath.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;

TAManyBodySDManager *TAManyBodySDManager::kInstance = nullptr;

TAManyBodySDManager::TAManyBodySDManager()
  : fManyBodySDListM(nullptr){}

TAManyBodySDManager *TAManyBodySDManager::Instance(){
  if(!kInstance) kInstance = new TAManyBodySDManager();
  return kInstance;
}

TAManyBodySDManager::~TAManyBodySDManager(){
  for(TAManyBodySD *&p : fManyBodySDVec){
    if(p){
      delete p; p = nullptr;
    } // end if
  } // end for
  fManyBodySDVec.clear();
} // end of the destructor

void TAManyBodySDManager::GenerateManyBodySD(){
  if(fManyBodySDVec.size()) return; // called already

  // obtain user input //
  string SPStatefile;
  short nParticle;
//  cout << "Please enter single particle state input file: ";
//  cin >> SPStatefile;
  SPStatefile = "sp.txt"; // DEBUG
//  cout << "Please enter number of particles: ";
  nParticle = 2; // DEBUG
//  if(!(cin >> nParticle)){
//    TAException::Error("TAManyBodySDManager",
//      "GenerateManyBodySD: Please enter an integer.");
//  }

  // generate SP state, MB state, and M-scheme MB state list //
  TASingleParticleStateManager *spStateManager
    = TASingleParticleStateManager::Instance();
  spStateManager->LoadSPListFile(SPStatefile);
  const int nSPState = spStateManager->GetSPStateVec().size();
  if(nParticle > nSPState){
    TAException::Error("TASingleParticleStateManager",
      "GenerateManyBodySD: The number of particles are larger than\
 the number of single particle states.");
  }
  const int nManyBodySD = TAMath::Binomial(nSPState, nParticle);
  fManyBodySDVec.clear(); fManyBodySDVec.reserve(nManyBodySD);


  /////////// odometer method to generate many-body basis /////////////
  int *SPStateVec = new int[nParticle];
  int index = 0; // many-body SD index
  // the first MBSD configuration
  for(int i = 0; i < nParticle; i++) SPStateVec[i] = i;
  fManyBodySDVec.push_back(new TAManyBodySD(index++, nParticle, SPStateVec));

  // generate the MBSDs //
  while(1){
    int i = nParticle - 1;
    int maxMBSD = i + nSPState - nParticle; // max MBSD of the i-th bit
    // get the bit that is not full
    while(i >= 0 && SPStateVec[i--] == maxMBSD--){}
    SPStateVec[++i]++; // increment the found unfull bit
    // recover the following bits
    while(i < nParticle - 1){
      SPStateVec[i + 1] = SPStateVec[i] + 1; i++;
    }
    fManyBodySDVec.push_back(new TAManyBodySD(index++, nParticle, SPStateVec));
    if(SPStateVec[0] == nSPState - nParticle) break;
  } // end while
  delete [] SPStateVec;
  ////////////////// END of the odometer algorithm /////////////////////

  if(!fManyBodySDVec.size())
    TAException::Error("TAManyBodySDManager",
      "GenerateManyBodySD: After called, still no ManyBodySD is generated.");
  if(int(fManyBodySDVec.size()) != TAMath::Binomial(nSPState, nParticle)){
    TAException::Error("TAManyBodySDManager",
      "GenerateManyBodySD: After called, number of ManyBodySD is not right.");
  }

  // display the geneated many-body SD for debugging purposes
  TAManyBodySD::kNParticle = nParticle;
  TAException::Info("TAManyBodySDManager",
    "GenerateManyBodySD: \n\
Display the generated many-body Slater determinants ~");
  for(TAManyBodySD *mp : fManyBodySDVec) mp->Print(); // DEBUG
  for(TAManyBodySD *mp : fManyBodySDVec) mp->PrintInBit(); // DEBUG
  cout << "Totally there're " << fManyBodySDVec.size();
  cout << " many-body Slater determinants in the list." << endl;
} // end member function GenerateManyBodySD

// generate the M-scheme many-body state basis
void TAManyBodySDManager::MSchemeGo(){
  if(fManyBodySDListM) return; // alrady called

  GenerateManyBodySD(); // Generate all the many-body basis

  // get the boundary of the MBSDs' M
  short min2M = (*std::min_element(fManyBodySDVec.begin(), fManyBodySDVec.end(),
    [](TAManyBodySD *a, TAManyBodySD *b){ return a->Get2M() < b->Get2M(); }))->Get2M();
  short max2M = (*std::max_element(fManyBodySDVec.begin(), fManyBodySDVec.end(),
    [](TAManyBodySD *a, TAManyBodySD *b){ return a->Get2M() < b->Get2M(); }))->Get2M();

  short twoM = -999;
//  cout << "Please enter 2M (total M*2 for the many-body system): ";
//  cin >> twoM;
  twoM = 0; // DEBUG
  if(twoM < min2M || twoM > max2M)
    TAException::Error("TAManyBodySDManager",
      "MSchemeGo: Input 2M: %d is not within [%d, %d]", twoM, min2M, max2M);

  // select the MBSDs with the M value specified by users //
  fManyBodySDListM = new TAManyBodySDList(twoM);
  for(TAManyBodySD *p : fManyBodySDVec){
    if(twoM == p->Get2M()) fManyBodySDListM->Add(p);
  }
  if(fManyBodySDListM->GetNBasis() == 0){
    TAException::Warn("TAManyBodySDManager",
      "MschemeGo: fManyBodySDListM is empty in the end.");
  }
  fManyBodySDListM->Print(); // DEBUG
  fManyBodySDListM->PrintInBit(); // DEBUG
} // end of member function MSchemeGo

TAManyBodySDList *TAManyBodySDManager::GetMBSDListM(){
  if(!fManyBodySDListM || !fManyBodySDVec.size()) MSchemeGo();
  return fManyBodySDListM;
}
