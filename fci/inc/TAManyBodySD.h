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

class TASPState;
class TAManyBodySDManager;

class TAManyBodySD{
public:
	TAManyBodySD(unsigned long long bit); ///< fPhase is set to 1
	virtual ~TAManyBodySD();

	short Get2M() const; ///< \retval the total jz*2
	double GetEnergy() const; ///< \retval the total sp energy
	const TABit &Bit() const{ return fBit; }

	void Print() const; ///< self-display
	void PrintInBit() const; ///< Print the many-body state in bit mode

	void GetSPStateArr(int *p) const; ///< return the spstae array
	TASPState *operator[](int i) const;
	int GetNParticle() const{ return fBit.count(); }
	/// note that this method works only if piared states are next to each other,
	/// and SPStates are ordered in ManyBodySD
	bool IsPaired() const; ///< \retval if there're broken pairs or not
	void SetBit(unsigned long long bit){ fBit = bit; }

	friend class TAManyBodySDManager;

protected:
	TABit fBit; ///< bit representation of this
};

#endif
