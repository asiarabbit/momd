/**
  MOMD project, Anyang Normal University, IMP-CAS
  \file TAMath.h
  \class TAMath
  \brief Math class, to provide some general math methods. Note that this is a
  tool class, so it is defined to be a static class.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/07/09
  \date Last modified: 2020/07/09 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright MOMD project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMath_h_
#define _TAMath_h_

class TAMath{
public:
  ~TAMath(){}

  static constexpr double Pi(){ return 3.14159265358979323846; }
  /// golden cut ratio
	static constexpr double Alpha(){ return 0.61803398874989484820458683436565; }
  /// rad per degree
	static constexpr double DEGREE(){ return 0.01745329251994329547; }
	static constexpr double Sqrt3(){ return 1.73205080756887719318; }

  /// sign of a number
  static double sign(double c);
};

#endif
