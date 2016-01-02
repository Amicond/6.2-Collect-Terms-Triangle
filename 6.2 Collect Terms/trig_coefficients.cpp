#include "stdafx.h"
#include "trig_coefficients.h"

const std::string trig_coefficients::angleName = "be";

trig_coefficients::	trig_coefficients(std::string angleName) //set coefficients name
{
	for (int i = 0; i < PowerAmount; i++)
		trig_angle_power.push_back(0);
}
bool trig_coefficients:: operator==(const trig_coefficients tc2)const
{
	for (unsigned int i = 0; i < trig_angle_power.size(); i++)
		if (trig_angle_power[i] != tc2.trig_angle_power[i])
			return false;
	return true;
}
bool trig_coefficients:: operator<(const trig_coefficients tc2)const
{
	for (unsigned int i = 0; i < PowerAmount; i++)
	{
		if (trig_angle_power[i] < tc2.trig_angle_power[i])
			return true;
		else if (tc2.trig_angle_power[i] < trig_angle_power[i])
			return false;
	}
	return false;
}
void trig_coefficients:: printTrigCoeff(std::ofstream &F)const
{
	switch (trig_angle_power[0]){
	case 0: break;
	case 1: F << "*Cos[" + angleName + "]"; break;
	default: F << "*Cos[" + angleName + "]^" << trig_angle_power[0];
	}
	switch (trig_angle_power[1]) {
	case 0: break;
	case 1:F << "*Sin[" + angleName + "]"; break;
	default:F << "*Sin[" + angleName + "]^" << trig_angle_power[1];
	}
	switch (trig_angle_power[2]){
	case 0: break;
	case 1: F << "*Cos[" + angleName + "/2]"; break;
	default: F << "*Cos[" + angleName + "/2]^" << trig_angle_power[2];
	}
	switch (trig_angle_power[3]){
	case 0: break;
	case 1:F << "*Sin[" + angleName + "/2]"; break;
	default:F << "*Sin[" + angleName + "/2]^" << trig_angle_power[3];
	}
}
void trig_coefficients::incCoeff(int type, int val)
{
	trig_angle_power[type] += val;
}
void trig_coefficients::clear()
{
	for (int i = 0; i < PowerAmount; i++)
		trig_angle_power[i] = 0;
}
std::string trig_coefficients::getPowers() const//for hash genertion
{
	std::string s = " ";
	for (auto elem : trig_angle_power)
		s += elem + " ";
	return s;
}
double trig_coefficients::return_coeff(double angle)
{
	double coeff = 1;

	for (int i = 0; i < trig_angle_power[0]; i++)
		coeff *= cos(angle);
	for (int i = 0; i < trig_angle_power[1]; i++)
		coeff *= sin(angle);
	for (int i = 0; i < trig_angle_power[2]; i++)
		coeff *= cos(angle / 2);
	for (int i = 0; i < trig_angle_power[3]; i++)
		coeff *= sin(angle / 2);
	return coeff;
}
