/**
	\file Hamiltonian.cxx
	\brief This code is for tackling the full configuration interaction
	in many-body problems which is a common question in many-body problems
*/

#include <iostream>
#include "TAHPairing.h"
#include "TADiagonalize.h"

using std::cout;
using std::endl;

int main(){
	TAHPairing *pair = new TAHPairing("../config/input.txt");
	matrix *mp = pair->Matrix();
	pair->PrintMBSD();
	mp->Print();

	const int n = pair->GetNMBSD(); // number of Slater determinants
	double d[n]; // , x[n]{}
	// x[0] = 1.;
	// for(auto &t : x) t = 1.;
	matrix v(n, n);
	TADiagonalize::JacobiSort(*mp, n, d, v);
	// TADiagonalize::TridiagQLSort(*mp, n, d);
	// TADiagonalize::JacobiSort(*mp, n, d);
	// TADiagonalize::Lanczos(*mp, n, d, v);
	for(int i = 0; i < n; i++) std::cout << d[i] << std::endl;
	mp->Print();

	return 0;
} // end of the main function
