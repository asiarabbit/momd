/**
	\file Hamiltonian.cxx
	\brief This code is for tackling the full configuration interaction
	in many-body problems which is a common question in many-body problems
*/

#include <iostream>
#include "TAHPairing.h"
#include "TADiagonalize.h"
#include "TATreeCol.h"

using std::cout;
using std::endl;

int main(){
	TAHPairing *pair = new TAHPairing("../config/input.txt");
	const matrix *H = pair->Matrix();
	H->Print();

	TATreeCol::CloseFile();
	return 0;
} // end of the main function
