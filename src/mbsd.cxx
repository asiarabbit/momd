/**
	SUNNY Project, Anyang Normal University, IMP-CAS
	\file mbsd.cxx
	\brief This code is for tackling the full configuration interaction
	in many-body problems which is a common question in many-body problems
	\date Created: 2020/01/31
	\date Last modified: 2020/01/31
	\copyright 2020, SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <iostream>
#include "TAManyBodySDManager.h"
#include "TADiagonalize.h"
#include "TAHPairing.h"
#include "TATreeCol.h"

using std::endl;
using std::cout;

int main(){
	// TAManyBodySDManager *mbsdManager = TAManyBodySDManager::Instance();
	// mbsdManager->LoadConfigFile("../config/input.txt");
	// mbsdManager->MSchemeGo(); // Generate the M-scheme many-body basis

	TAHPairing *pair = new TAHPairing("../config/input.txt");
	matrix *H = pair->Matrix();
	const int n = pair->GetNBasis();

	matrix z(n, n);
	double d[n];
	// TADiagonalize::LanczosPurge(*H, n, d, z);
	TADiagonalize::JacobiSort(*H, n, d, z);
	TADiagonalize::TridiagQLSort(*H, n, d);

	for(int i = 0; i < n; i++) cout << d[i] << endl;
	H->Print();

	TATreeCol::CloseFile();
	return 0;
}
