/**
	\file Hamiltonian.cxx
	\brief This code is for tackling the full configuration interaction
	in many-body problems which is a common question in many-body problems
*/

#include "TAFCI.h"

int main(){
	TAFCI *fci = TAFCI::Instance();
	fci->Go(); // Generate M-scheme many-body basis and establish H matrix

	return 0;
} // end of the main function
