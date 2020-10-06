/**
	\file TASingleParticleState.h
	\class TASingleParticleState
	\brief Single particle state for a one-body Hamiltonian with a central force.
	Labeled by n, l, j and mj. Used as a building basis for configuration
	interaction calculations in quantum many-body problems.
	\author SUN Yazhou, asia.rabbit@163.com
	\date Created: 2020/01/31
	\date Last modified: 2020/01/31 by SUN Yazhou
	\copyright 2020, SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TASingleParticleState_h_
#define _TASingleParticleState_h_

class TASingleParticleState{
public:
	TASingleParticleState(int index, short n, short l,
		short _2j, short _2mj, double energy);
	virtual ~TASingleParticleState();
	short GetMj(){ return f2mj; } /// \return the third component mj
	double GetEnergy(){ return fEnergy; }
	void Print() const; ///< print the single particle state

protected:
	int fIndex; ///< index of the single particle state
	short fn; ///< number of nodes (zero-points) for the radial functiion
	short fl; ///< orbital angular momentum number
	short f2j; ///< angular momentum number*2
	short f2mj; ///< the third compoment jz*2
	double fEnergy; ///< the energy of the single-particle state
};

#endif
