/**
	SUNNY Project, Anyang Normal University, IMP-CAS
	\file TAManyBodySDList.h
	\class TAManyBodySDList
	\brief A list to store many-body Slater determinants, classified by M, which is
	the 3rd component jz of the total angular momentum. So M is the same for each
	member of this list.
	\author SUN Yazhou
	\date Created: 2020/01/31
	\date Last modified: 2020/02/11 by SUN Yazhou
	\copyright 2020 SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAManyBodySDList_h_
#define _TAManyBodySDList_h_

#include "TAMBSDTree.h"

class TASPState;

using std::vector;

class TAManyBodySDList{
public:
	TAManyBodySDList(short twoM);
	virtual ~TAManyBodySDList();

	short Get2M() const{ return f2M; }
	TAMBSDTree *GetMBSDTree(){ return fManyBodySDTree; }
	bool MBSDFilter(unsigned long long bit) const; ///< filter mbsds
	void Add(unsigned long long bit);
	void Add(int *spsArr, int nParticle);
	void Print();
	void PrintInBit(); ///< Print all the mbsd-s in bit mode
	unsigned long long GetNBasis() const;
	unsigned long long Bit(unsigned long long mbsdIndex);
	unsigned long long operator[](unsigned long long mbsdIndex){ return Bit(mbsdIndex); }
	TASPState *GetSPState(unsigned long long mbsdIndex, int i);
	void Save(){ fManyBodySDTree->Save(); } ///< save the many-body state tree

	/// \retval <rr|a+_p * a_q|cc>, 1N force
	int Integral(int rr, int p, int q, int cc);
	/// \retval <rr|a+_p*a+_q * a_s*a_r|cc>, 2N force
	int Integral(int rr, int p, int q, int s, int r, int cc);
	/// \retval <rr|a+_p*a+_q*a+_r * a_u*a_t*a_s|cc>, 3N force
	int Integral(int rr, int p, int q, int r, int u, int t, int s, int cc);

protected:
	short f2M; ///< the uniform M*2 for this list
	TAMBSDTree *fManyBodySDTree;
};

#endif
