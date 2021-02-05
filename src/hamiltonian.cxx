/**
	\file Hamiltonian.cxx
	\brief This code is for tackling the full configuration interaction
	in many-body problems which is a common question in many-body problems
*/

#include "TAHPairing.h"

int main(){
	TAHPairing *pair = new TAHPairing("../config/input.txt");
	matrix *mp = pair->Matrix();
	mp->Print();


	return 0;
} // end of the main function
