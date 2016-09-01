#include "stdafx.h"
#include "ground_energy.h"
#include "cos.h" 



groundEnergy::trigPowers::trigPowers() {
	for (int i = 0; i < Types; i++)
	{
		sinPower[i] = 0;
		cosPower[i] = 0;
	}
}

void groundEnergy::trigPowers::add_trig(int type, int power, int plaquet_type)
{
	switch (type)
	{
		case 0: sinPower[plaquet_type]++; break;
		case 1: cosPower[plaquet_type]++; break;
		default: std::cout << "\n\nTRIG POWER ALARM!\n\n"; std::cin >> type;
	}
}

bool groundEnergy::trigPowers:: operator==(const trigPowers& right)const
{
	bool res = true;
	for (int i = 0; i < Types; i++)
	{
		if (sinPower[i] != right.sinPower[i]||cosPower[i]!=right.cosPower[i])
		{
			res = false;
			break;
		}
	}
	return res;
}

void groundEnergy::set(double factor, bool IfNumerical)
{
	this->spin = factor;
	this->ifNumerical = IfNumerical;
}
void groundEnergy::addTerm(unsigned int order, term t1)
{
	while (energy.size() <= order)
		energy.push_back(0);
	if (t1.len == -1)// C-case
	{
		energy[order] += t1.value;
	}
	else
	{
		double res = t1.value;
		for (unsigned int i = 0; i < t1.len; i++)
			res *= spin;
		res /= t1.len;
		energy[order] += res;
	}
}
void groundEnergy::clearTerms()
{
	energy.clear();
	for (int i = 0; i < N; i++)
		trigEnergy[i].clear();
}

double groundEnergy::returnE(unsigned int order)
{
	if ( order <= energy.size())
	{
		return energy[order];
	}
	else
		return 99999;
}
//rotation case


//void groundEnergy::addTermRotationZeroPi2Angles(int order, term t1)
//{
//	bool flag = false;//no such term by efault
//	trigPowers cur;
//	cur.value = t1.value;
//	bool second_type;
//	if (t1.len != -1)// not C-case
//	{
//		cur.value /= t1.len;
//		for (unsigned int i = 0; i < t1.len; i++){
//			second_type = Cos::getSignZeroPi(t1.nums[i]) == -1;
//			if (t1.ops[i] == 'z')
//				cur.add_trig(0, 1, second_type);
//			else 
//				cur.add_trig(1, 1, second_type);
//		}
//	}
//
//	for (unsigned int i = 0; i < trigEnergy[order].size(); i++)
//	{
//		if (trigEnergy[order][i] == cur)
//		{
//			trigEnergy[order][i].value += cur.value;
//			flag = true;
//			break;
//		}
//	}
//	if (!flag)
//		trigEnergy[order].push_back(cur);
//}
//
//void groundEnergy::addTermRotationPiZero2Angles(int order, term t1)
//{
//	bool flag = false;//no such term by efault
//	trigPowers cur;
//	cur.value = t1.value;
//	bool second_type;
//	if (t1.len != -1)// not C-case
//	{
//		cur.value /= t1.len;
//		for (unsigned int i = 0; i < t1.len; i++){
//			second_type = Cos::getSignPiZero(t1.nums[i]) == -1;
//			if (t1.ops[i] == 'z')
//				cur.add_trig(0, 1, second_type);
//			else
//				cur.add_trig(1, 1, second_type);
//		}
//	}
//
//	for (unsigned int i = 0; i < trigEnergy[order].size(); i++)
//	{
//		if (trigEnergy[order][i] == cur)
//		{
//			trigEnergy[order][i].value += cur.value;
//			flag = true;
//			break;
//		}
//	}
//	if (!flag)
//		trigEnergy[order].push_back(cur);
//}
//
//void groundEnergy::addTermRotation2sublattices(int order, term t1)
//{
//	bool flag = false;//no such term by efault
//	trigPowers cur;
//	cur.value = t1.value;
//	cur.sinPower = 0;
//	cur.cosPower = 0;
//
//	if (t1.len != -1)// not C-case
//	{
//		cur.value /= t1.len;
//		for (unsigned int i = 0; i < t1.len; i++)
//			if (t1.ops[i] == 'z'){
//				cur.cosPower++;
//			}
//			else {
//				cur.sinPower++;
//				cur.value *= Cos::getSign(t1.nums[i]);
//			}
//	}
//
//	for (unsigned int i = 0; i < trigEnergy[order].size(); i++)
//	{
//		if (trigEnergy[order][i] == cur)
//		{
//			trigEnergy[order][i].value += cur.value;
//			flag = true;
//			break;
//		}
//	}
//	if (!flag)
//		trigEnergy[order].push_back(cur);
//}


void groundEnergy::addTriangleTermRotation(int order, term t1,int zeroType)
{
	bool flag = false;//no such term by efault
	trigPowers cur;
	cur.value = t1.value;

	if (t1.len != -1)// not C-case
	{
		for (unsigned int i = 0; i < t1.len; i++) 
		{
			if (t1.ops[i] == 'z')
				cur.add_trig(1, 1,pseudospin_types::getType(t1.nums[i], zeroType));
			else
				cur.add_trig(0, 1, pseudospin_types::getType(t1.nums[i], zeroType));
		}
	}

	for (unsigned int i = 0; i < trigEnergy[order].size(); i++)
	{
		if (trigEnergy[order][i] == cur)
		{
			trigEnergy[order][i].value += cur.value;
			flag = true;
			break;
		}
	}
	if (!flag)
		trigEnergy[order].push_back(cur);
}


void groundEnergy::printTermRotation(std::ostream& out, int order,bool two_angels)
{
	for (unsigned int i = 0; i < trigEnergy[order].size(); i++)
	{
		
		out << trigEnergy[order][i].value;
		int sum_power=0;
		for (int j = 0; j < Types; j++) {

			switch (trigEnergy[order][i].cosPower[j]) {
			case 0: break;
			case 1: out << "*Cos["<< pseudospin_types::getAngleName(j) <<"]"; break;
			default:out << "*Cos[" << pseudospin_types::getAngleName(j)<<"]^" << trigEnergy[order][i].cosPower[j]; break;
			}
			switch (trigEnergy[order][i].sinPower[j]) {
			case 0: break;
			case 1: out << "*Sin[" << pseudospin_types::getAngleName(j) << "]"; break;
			default: out << "*Sin[" << pseudospin_types::getAngleName(j) << "]^" << trigEnergy[order][i].sinPower[j]; break;
			}
			sum_power += trigEnergy[order][i].cosPower[j] + trigEnergy[order][i].sinPower[j];
		};
		

		if (!ifNumerical){
			switch (sum_power){
				case 0: break;
				case 1: out << "*Sz"; break;
				default: out << "*Sz^" << (sum_power);
			}
		}
		else {
			switch (sum_power){
				case 0: break;
				case 1: out << "*(" << spin << ")"; break;
				default: out << "*(" << spin << ")^" << (sum_power);
			}
			
		}
	
		if (i != trigEnergy[order].size() - 1)
			out << "+";
	}
}
