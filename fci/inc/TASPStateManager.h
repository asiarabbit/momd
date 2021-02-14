/**
	\file TASPStateManager.h
	\class TASPStateManager
	\brief A list of TASingleParticle objects, responsible for reading the user-input
	single-particle (SP) state files, generating SP state objects, and managing the SP
	states. Note that this is a singleton class.
	\author SUN Yazhou
	\date Created: 2020/01/31
	\date Last modified: 2020/01/31 by SUN Yazhou
	\copyright 2020 SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASPStateManager_h_
#define _TASPStateManager_h_

#include <vector>
#include <string>

using std::vector;
using std::string;

class TASPStateManager{
public:
	virtual ~TASPStateManager();
	static TASPStateManager *Instance();
	vector<TASPState *> &GetSPStateVec(){ return fSPStateVec; }
	int GetNSPState() const;
	/// \param file: the input file is of format as follows:
	/// index n l 2j 2mj energy
	/// Lines starting with # are ignored
	void LoadSPListFile(const string &file);

protected:
	TASPStateManager();
	static TASPStateManager *kInstance;
	string fFileIn;
	vector<TASPState *> fSPStateVec;
};

#endif
