/**
	MOMD project, Anyang Normal University, IMP-CAS
	\file TAException.cxx
	\class TAException
	\brief to print error, warnings, etc.
	\author SUN Yazhou, asia.rabbit@163.com
	\date Created: 2020/07/08
	\date Last modified: 2020/10/08
	\copyright 2020, by SUN Yazhou
	\copyright MOMD project, Anyang Normal University, IMP-CAS
*/


#include <iostream>
#include <ctime>
#include <stdarg.h> // variable arguments handle function declarations
#include <cstdlib>
#include "TAException.h"

using std::cout;
using std::endl;

// assignment and definition
bool TAException::fIsSilent = false;
bool TAException::fIsDebug = true;

void TAException::Info(const char *cname, const char *fmt, ...){
	if(fIsSilent) return;

	va_list arg_ptr;
	va_start(arg_ptr, fmt);
	char msg[1024], omsg[1024];
	sprintf(msg, "  SUNNY-Info: \033[32;1m[%s]\033[33;1m %s\n\033[0m",
		cname, fmt);
	vsprintf(omsg, msg, arg_ptr);
	va_end(arg_ptr);
	cout << omsg;
} // end of member function Info

void TAException::Error(const char *cname, const char *fmt, ...){
	static bool nomore = false;
	if(nomore) return;

	va_list arg_ptr;
	va_start(arg_ptr, fmt);
	char msg[1024], omsg[1024];
	sprintf(msg, "\033[31;1m  SUNNY-Error: \033[0m\033[33;1m[%s]\033[0m\
\033[31;1m %s\n\033[0m", cname, fmt);
	vsprintf(omsg, msg, arg_ptr);
	va_end(arg_ptr);
	cout << omsg;

	// select repeatability
	if(!nomore){
		cout << "[q]: abort, [n]: do not stop any more, others: continue > "
<< std::flush;
		std::string str;
		std::getline(std::cin, str);
		if(!str.empty() && (str.c_str())[0] == 'q'){
			cout << "Aborting SUNNY...\n";
			exit(1);
		}
		else if(!str.empty() && (str.c_str())[0] == 'n') nomore = true;
	}
} // end of member function Error

void TAException::Warn(const char *cname, const char *fmt, ...){
	static bool nomore = false;
	if(nomore) return;

	va_list arg_ptr;
	va_start(arg_ptr, fmt);
	char msg[1024], omsg[1024];
	sprintf(msg, "\033[36;1m  SUNNY-Warning: \033[0m\033[33;1m[%s]\
\033[0m\033[36;1m %s\n\033[0m", cname, fmt);
	vsprintf(omsg, msg, arg_ptr);
	va_end(arg_ptr);
	cout << omsg;

	// select repeatability
	if(!nomore){
		cout << "[q]: abort, [n]: do not stop any more, others: continue > "
<< std::flush;
		std::string str;
		std::getline(std::cin, str);
		if(!str.empty() && (str.c_str())[0] == 'q'){
			cout << "Aborting SUNNY...\n";
			exit(1);
		}
		else if(!str.empty() && (str.c_str())[0] == 'n') nomore = true;
	} // end if
} // end of member function Warn

void TAException::Debug(const char *cname, const char *fmt, ...){
	if(!fIsDebug) return;

	va_list arg_ptr;
	va_start(arg_ptr, fmt);
	char msg[1024], omsg[1024];
	sprintf(msg, "  SUNNY-Debug: \033[33;1m[%s]\033[0m %s\n", cname, fmt);
	vsprintf(omsg, msg, arg_ptr);
	va_end(arg_ptr);
	cout << omsg;

	cout << "Press ENTER to continue...";
	getchar();
} // end of member function Debug
