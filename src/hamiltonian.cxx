/**
	\file Hamiltonian.cxx
	\brief This code is for tackling the full configuration interaction
	in many-body problems which is a common question in many-body problems
*/

#include <iostream>
#include "TAHPairing.h"
#include "TADiagonalize.h"
#include "TALanczos.h"
#include "TASparseMatrix.h"

using std::cout;
using std::endl;

int main(){
	TAHPairing *pair = new TAHPairing("../config/input.txt");
	pair->Matrix()->Print();

	const short n = pair->GetNBasis() < 30 ? pair->GetNBasis() : 30;
	TASparseMatrix z; double d[n*2]{};
	TALanczos::PurgeRestart(*pair, d, z, n);

	for(int i = 0; i < n; i++) cout << d[i] << endl;
	z.Print();

	z.Save();
	TATreeCol::CloseFile();
	return 0;
} // end of the main function
