/**
	\file TASPState.C
	\class TASPState
	\brief Single particle state for a one-body Hamiltonian with a central force.
	Labeled by n, l, j and mj. Used as a building basis for configuration
	interaction calculations in quantum many-body problems.
	\author SUN Yazhou, asia.rabbit@163.com
	\date Created: 2020/01/31
	\date Last modified: 2020/01/31 by SUN Yazhou
	\copyright 2020, SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include <iomanip>
#include "TASPState.h"

using std::cout;
using std::setw;

TASPState::	TASPState(int index, short n, short l, short _2j, short _2mj, double energy)
		: fIndex(index), fn(n), fl(l), f2j(_2j), f2mj(_2mj), fEnergy(energy){}

TASPState::~TASPState(){}

/// print the single particle state
void TASPState::Print() const{
	cout << std::right;
	cout << "index: " << setw(3) << fIndex;
	cout << "   n: " << setw(2) << fn << "   l: " << setw(1) << fl;
	cout << "   2j: " << setw(1) << f2j << "   2mj: " << setw(3) << f2mj;
	cout << "   energy: " << setw(3) << fEnergy << std::endl;
	cout << std::left;
} // end of member function Print
