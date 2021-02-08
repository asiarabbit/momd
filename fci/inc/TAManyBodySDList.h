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

#include <vector>

class TAManyBodySD;

using std::vector;

class TAManyBodySDList{
public:
	TAManyBodySDList(short twoM);
	virtual ~TAManyBodySDList();

	short Get2M() const{ return f2M; }
	vector<TAManyBodySD *> &GetManyBodySDVec(){ return fManyBodySDVec; }
	void Add(TAManyBodySD *mbsd);
	void Pairing(); /// remove those with broken pairs, i.e. >2 single particles
	void Print() const;
	void PrintInBit() const; ///< Print all the mbsd-s in bit mode
	int GetNBasis() const{ return fManyBodySDVec.size(); }
	TAManyBodySD *operator[](int i) const; ///< \retval fManyBodySDVec[i]

	/// \retval <rr|a+_p * a_q|cc>, 1N force
	int Integral(int rr, int p, int q, int cc) const;
	/// \retval <rr|a+_p*a+_q * a_s*a_r|cc>, 2N force
	int Integral(int rr, int p, int q, int s, int r, int cc) const;
	/// \retval <rr|a+_p*a+_q*a+_r * a_u*a_t*a_s|cc>, 3N force
	int Integral(int rr, int p, int q, int r, int u, int t, int s, int cc) const;

protected:
	short f2M; ///< the uniform M*2 for this list
	vector<TAManyBodySD *> fManyBodySDVec;
};

#endif
