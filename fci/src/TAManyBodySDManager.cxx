/**
  \file TAManyBodySDManager.cxx
  \class TAManyBodySDManager
  \brief A class to generate many-body basis and manage TAManyBOdySDList objects.
  \author SUN Yazhou
  \date Created: 2020/02/01
  \date Last modified: 2020/02/02 by SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include "TAManyBodySDManager.h"
#include "TAManyBodySDList.h"
#include "TAManyBodySD.h"
#include "TASPStateManager.h"
#include "TAException.h"

using std::cout;
using std::cin;
using std::endl;

TAManyBodySDManager *TAManyBodySDManager::kInstance = nullptr;

TAManyBodySDManager::TAManyBodySDManager() : fManyBodySDListM(nullptr),
    fNParticle(0), fNSPState(0), f2M(-9999){
  fSPStatefile = "";
} // end of the constructor

TAManyBodySDManager *TAManyBodySDManager::Instance(){
  if(!kInstance) kInstance = new TAManyBodySDManager();
  return kInstance;
}

TAManyBodySDManager::~TAManyBodySDManager(){}

// skip spaces and tabs, return subscript of the valid char
inline int skipCrap(const char *s){
	int tmp = 0;
	while(1){
		char c = s[tmp++];
		if(' ' != c && '\t' != c) break;
	}
	return tmp - 1;
}
/// assign nparticle, nspstate, 2*M and the spstate space
void TAManyBodySDManager::LoadConfigFile(const string &file){
  std::ifstream inFile(file);
	if(!inFile.is_open()){
		TAException::Error("TAManyBodySDManager",
			"LoadConfigFile: input file %s open error.", file);
	}
	char line[512] = "";
  std::stringstream ss;
	while(inFile.getline(line, sizeof(line))){
		if(0 == strlen(line)) continue; // blank line
		int tmp = skipCrap(line);
		if('#' == line[tmp] || '\0' == line[tmp]) continue; // commentary line
    ss << line; ss << ' ';
	} // end while
  ss >> fSPStatefile; ss >> fNParticle; ss >> f2M;

	inFile.close();
  cout << "Single particle state input file: " << endl;
  cout << fSPStatefile << endl;
  cout << "Number of particles: " << fNParticle << endl;
  cout << "2M for M-scheme: " << f2M << endl << endl;
} // end of member function LoadConfigFile


void TAManyBodySDManager::GenerateManyBodySD(){
  if(fManyBodySDListM->GetNBasis()) return; // called already

  // generate SP state, MB state, and M-scheme MB state list //
  if("" == fSPStatefile) TAException::Error("TAManyBodySDManager",
    "GenerateManyBodySD: Single-particle state space inputfile not assigned yet.");
  TASPStateManager *spStateManager = TASPStateManager::Instance();
  spStateManager->LoadSPListFile(fSPStatefile);
  fNSPState = spStateManager->GetNSPState();
  if(fNParticle > fNSPState) TAException::Error("TAManyBodySDManager",
    "GenerateManyBodySD: nParticle > nSPState.");

  /////////// odometer method to generate many-body basis /////////////
  if(0 == fNParticle) TAException::Error("TAManyBodySDManager",
    "GenerateManyBodySD: Number of particles not assigned.");
  int *SPStateVec = new int[fNParticle];
  // the first MBSD configuration
  for(int i = 0; i < fNParticle; i++) SPStateVec[i] = i;
  fManyBodySDListM->Add(SPStateVec, fNParticle); ///< Fill mbsd in a TTree obj

  // generate the many-body slater determinants //
  // fNParticle sp-states are stored in each int of array SPStateVec. The sp-states are
  // labelled from 0, and all the states are distinct, e.g. [0,1,2,3], [0,1,2,4], ...
  // the odometer method starts from the last sp-state, increments it to the maximum
  // (maxSP), and then recover the following sp-states. e.g. [3,4,5,6]->[4,5,6,7]
  // for 4 particles in 8 sp-states.
  while(1){
    int i = fNParticle - 1, maxSP = fNSPState - 1; // max state of the last particle
    // find the bit that is not full (overflow) in array SPStateVec backwards
    while(i >= 0 && SPStateVec[i--] == maxSP--);
    SPStateVec[++i]++; // increment the found unfull bit
    // recover the following (rightward) sp-states
    while(i < fNParticle - 1){
      SPStateVec[i + 1] = SPStateVec[i] + 1; i++;
    }
    fManyBodySDListM->Add(SPStateVec, fNParticle); ///< Fill mbsd in a TTree obj
    if(fNSPState - fNParticle == SPStateVec[0]) break;
  } // end while
  delete [] SPStateVec;
  ////////////////// END of the odometer algorithm /////////////////////

  unsigned long long nb = fManyBodySDListM->GetNBasis();
  if(!nb) TAException::Error("TAManyBodySDManager",
      "GenerateManyBodySD: After called, still no ManyBodySD is generated.");

  // display the geneated many-body SD for debugging purposes
  TAException::Info("TAManyBodySDManager",
    "GenerateManyBodySD: Display the generated many-body Slater determinants ~");
  cout << "Totally there're " << fManyBodySDListM->GetNBasis(); // DEBUG
  cout << " many-body Slater determinants in the list" << endl; // DEBUG
} // end member function GenerateManyBodySD

// generate the M-scheme many-body state basis
void TAManyBodySDManager::MSchemeGo(){
  if(fManyBodySDListM) return; // already called

  if(-9999 == f2M) TAException::Error("TAManyBodySDManager",
    "GenerateManyBodySD: f2M for M-scheme not assigned.");
  fManyBodySDListM = new TAManyBodySDList(f2M);
  GenerateManyBodySD(); // Generate all the many-body basis, plus filter

  if(0 == fManyBodySDListM->GetNBasis()) TAException::Warn("TAManyBodySDManager",
    "MSchemeGo: fManyBodySDListM is empty in the end.");
  // fManyBodySDListM->Print(); // DEBUG
  // fManyBodySDListM->PrintInBit(); // DEBUG
  fManyBodySDListM->Save(); // save the mbsd TTree obj
} // end of member function MSchemeGo

TAManyBodySDList *TAManyBodySDManager::GetMBSDListM(){
  if(!fManyBodySDListM) TAException::Warn("TAManyBodySDManager",
    "GetMBSDListM: fManyBodySDListM nullptr. MSchemeGo() not called yet?");
  return fManyBodySDListM;
}
