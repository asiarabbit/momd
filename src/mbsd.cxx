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

#include "TAManyBodySDManager.h"

int main(){
	TAManyBodySDManager *mbsdManager = TAManyBodySDManager::Instance();
	mbsdManager->MSchemeGo(); // Generate the M-scheme many-body basis

	return 0;
}
