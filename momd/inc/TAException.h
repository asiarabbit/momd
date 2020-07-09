/**
	MOMD project, Anyang Normal University, IMP-CAS
	\file TAException.h
	\class TAException
	\brief to print error, warnings, etc.
	\author SUN Yazhou, asia.rabbit@163.com
	\date Created: 2020/07/08
	\date Last modified: 2020/07/08
	\copyright 2020, by SUN Yazhou
	\copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAException_h_
#define _TAException_h_

class TAException{
public:
	TAException(){};
	virtual ~TAException(){};

	/// \param cname: class name
	static void Info(const char *cname, const char *fmt, ...);
	static void Error(const char *cname, const char *fmt, ...);
	static void Warn(const char *cname, const char *fmt, ...);
	static void Debug(const char *cname, const char *fmt, ...);

	static void Silent(bool opt = true){ fIsSilent = opt; }
	static bool IsSilent() { return fIsSilent; }
	static void SetDebug(bool opt = true){ fIsDebug = opt; }
	static bool IsDebug() { return fIsDebug; }
protected:
	static bool fIsSilent; ///< to switch Info() method
	static bool fIsDebug; ///< to switch Debug() method
};

#endif
