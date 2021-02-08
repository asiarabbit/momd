/**
	\file TASingleParticleStateManager.h
	\class TASingleParticleStateManager
	\brief A list of TASingleParticle objects, responsible for reading the user-input
	single-particle (SP) state files, generating SP state objects, and managing the SP
	states. Note that this is a singleton class.
	\author SUN Yazhou
	\date Created: 2020/01/31
	\date Last modified: 2020/01/31 by SUN Yazhou
	\copyright 2020 SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASingleParticleStateManager_h_
#define _TASingleParticleStateManager_h_

#include <vector>
#include <string>

using std::vector;
using std::string;

class TASingleParticleStateManager{
public:
	virtual ~TASingleParticleStateManager();
	static TASingleParticleStateManager *Instance();
	vector<TASingleParticleState *> &GetSPStateVec(){ return fSPStateVec; }
	int GetNSPState() const;
	/// \param file: the input file is of format as follows:
	/// index n l 2j 2mj energy
	/// Lines starting with # are ignored
	void LoadSPListFile(const string &file);

protected:
	TASingleParticleStateManager();
	static TASingleParticleStateManager *kInstance;
	string fFileIn;
	vector<TASingleParticleState *> fSPStateVec;
};

#endif
