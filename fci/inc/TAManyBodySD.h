/**
	SUNNY Project, Anyang Normal University, IMP-CAS
	\file TAManyBodySD.h
	\class TAManyBodySD
	\brief Slater determinant (SD) class for many-body problems. Each SD represents
	a configuration of the nucleons in the single-particle state.
	\author SUN Yazhou
	\date Created: 2020/01/31
	\date Last modified: 2020/02/11 by Sun Yazhou
	\copyright 2020 SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAManyBodySD_h_
#define _TAManyBodySD_h_

#include "TABit.h"

class TASingleParticleState;
class TAManyBodySDManager;

class TAManyBodySD{
public:
	TAManyBodySD(int index, int nParticle, int *SPState);
	virtual ~TAManyBodySD();
	short Get2M() const{ return f2M; } ///< \retval the total jz*2
	const TABit &Bit() const{ return fBit; }
	void Print() const; ///< self-display
	void PrintInBit() const; ///< Print the many-body state in bit mode
	/// note that this method works only if piared states are next to each other,
	/// and SPStates are ordered in ManyBodySD
	bool IsPaired() const; ///< \retval if there're broken pairs or not

	TASingleParticleState *operator[](int i);
	void SetIndex(int index){ fIndex = index; }
	void UpdateSPStateArr(); ///< update fSPStateArr to fBit
	void UpdateBit(); ///< update fBit to fSPStateArr
	int *IntArr() const; ///< \retval fSPStateArr, update with fBit if necessary
	int GetNParticle() const;

	friend class TAManyBodySDManager;

protected:
	int fNParticle;
	int fIndex; ///< index of the object
	short f2M; ///< the total jz*2
	double fEnergy; ///< the total energy of the SD
	TABit fBit; // bit representation of this
	/// Dynamical memory allocation. The length is the number of particles
	int *fSPStateArr;
};

#endif
