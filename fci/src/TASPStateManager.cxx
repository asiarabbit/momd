/**
	\file TASPStateManager.C
	\class TASPStateManager
	\brief A list of TASPState objects, responsible for reading the
	 user-input single-particle (SP) state files, generating SP state objects,
	and managing the SP states. Note that this is a singleton class.
	\author SUN Yazhou
	\date Created: 2020/01/31
	\date Last modified: 2020/01/31 by SUN Yazhou
	\copyright 2020 SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <string>
#include <cstring>
#include <fstream>
#include "TASPState.h"
#include "TASPStateManager.h"
#include "TAException.h"

using std::ifstream;
using std::string;

TASPStateManager *TASPStateManager::kInstance = nullptr;

TASPStateManager::TASPStateManager(){}

TASPStateManager *TASPStateManager::Instance(){
	if(!kInstance) kInstance = new TASPStateManager();
	return kInstance;
}

TASPStateManager::~TASPStateManager(){
	for(TASPState *&p : fSPStateVec){
		if(p){
			delete p; p = nullptr;
		} // end if
	} // end for
	fSPStateVec.clear();
} // end of the destructor

int TASPStateManager::GetNSPState() const{
	if(!fSPStateVec.size()){
		TAException::Warn("TASPStateManager",
			"GetNSPState: The fSPStateVec is empty.");
	}
	return fSPStateVec.size();
} // end of member function GetNSPState

// skip spaces and tabs, return subscript of the valid char
inline int skipCrap(const char *s){
	int tmp = 0;
	while(1){
		char c = s[tmp++];
		if(' ' != c && '\t' != c) break;
	}
	return tmp - 1;
}
/// \param file: the input file is of format as follows:
/// index n l 2j 2mj energy
/// Lines starting with # are ignored
void TASPStateManager::LoadSPListFile(const string &file){
	fFileIn = file;
	const char *filename = file.c_str();
	ifstream ff(filename);
	if(!ff.is_open()){
		TAException::Error("TASPStateManager",
			"LoadSPListFile: input file %s open error.", file);
	}
	char line[512] = "";
	while(ff.getline(line, sizeof(line))){
		if(0 == strlen(line)) continue; // blank line
		int tmp = skipCrap(line);
		if('#' == line[tmp] || '\0' == line[tmp]) continue; // commentary line

		int index, n, l, two_j, two_mj; double energy;
		sscanf(line, "%d %d %d %d %d %lg", &index, &n, &l, &two_j, &two_mj, &energy);
		fSPStateVec.push_back(
			new TASPState(index, n, l, two_j, two_mj, energy));
	} // end while
	ff.close();

	// TAException::Info("TASPStateManager",
	// 	"LoadSPListFile: \nDisplay the SP States read in ~");
	// for(TASPState *sp : fSPStateVec) sp->Print(); // DEBUG
} // end of member function LoadSPListFile()
