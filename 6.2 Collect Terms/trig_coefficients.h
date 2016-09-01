#include "stdafx.h"

#ifndef __TRIG_COEFFS_H__
#define __TRIG_COEFFS_H__
class trig_coefficients
{
	//terms 
	// 0 : cos(theta)
	// 1 : sin(theta)
	// 2 : cos^2(theta/2)
	// 3 : sin^2(theta/2)
	static const std::string angleName;
	static const int PowerAmount = 4;
	std::vector<int> trig_angle_power;
public: 
	trig_coefficients(std::string angleName = "beta"); //set coefficients name
	
	bool operator==(const trig_coefficients tc2)const;
	
	bool operator<(const trig_coefficients tc2)const;
	
	void printTrigCoeff(std::ofstream &F)const;
	
	void incCoeff(int type, int val);
	
	void clear();
	
	std::string getPowers() const;//for hash genertion
	
	double return_coeff(double angle);
	
};

#endif //__TRIG_COEFFS_H__