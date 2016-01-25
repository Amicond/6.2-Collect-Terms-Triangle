#include "stdafx.h"
#include "ground_energy.h"
#include "cos.h"



groundEnergy::trigPowers::trigPowers() {
	complex_power = 0;
	sinPower = 0;
	cosPower = 0;
}

bool groundEnergy::trigPowers:: operator==(const trigPowers& right)const
{
	return (cosPower == right.cosPower) && (sinPower == right.sinPower)&&(complex_power==right.complex_power);
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
void groundEnergy::addTermRotation(int order, term t1)
{
	bool flag = false;//no such term by efault
	trigPowers cur;
	cur.value = t1.value;
	cur.sinPower = 0;
	cur.cosPower = 0;

	if (t1.len != -1)// not C-case
	{
		cur.value /= t1.len;
		for (unsigned int i = 0; i < t1.len; i++)
			if (t1.ops[i] == 'z')
				cur.cosPower++;
			else
				cur.sinPower++;
	}

	for (unsigned int i = 0; i < trigEnergy[order].size(); i++)
	{
		if (trigEnergy[order][i] == cur)
		{
			trigEnergy[order][i].value += cur.value;
			flag = true;
		}
	}
	if (!flag)
		trigEnergy[order].push_back(cur);
}

void groundEnergy::addTermRotationAntiferromagnet(int order, term t1)
{
	bool flag = false;//no such term by efault
	trigPowers cur;
	cur.value = t1.value;
	cur.sinPower = 0;
	cur.cosPower = 0;

	if (t1.len != -1)// not C-case
	{
		cur.value /= t1.len;
		for (unsigned int i = 0; i < t1.len; i++)
			if (t1.ops[i] == 'z'){
				cur.cosPower++;
				cur.value *= Cos::getSign(t1.nums[i]);
			}
			else
				cur.sinPower++;
	}

	for (unsigned int i = 0; i < trigEnergy[order].size(); i++)
	{
		if (trigEnergy[order][i] == cur)
		{
			trigEnergy[order][i].value += cur.value;
			flag = true;
		}
	}
	if (!flag)
		trigEnergy[order].push_back(cur);

}

void groundEnergy::addTermRotationPiZero(int order, term t1)
{
	bool flag = false;//no such term by efault
	trigPowers cur;
	cur.value = t1.value;
	cur.sinPower = 0;
	cur.cosPower = 0;

	if (t1.len != -1)// not C-case
	{
		cur.value /= t1.len;
		for (unsigned int i = 0; i < t1.len; i++)
			if (t1.ops[i] == 'z'){
				cur.cosPower++;
			}
			else {
				cur.sinPower++;
				cur.value *= Cos::getSignPiZero(t1.nums[i]);
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

void groundEnergy::addTermRotation2sublattices(int order, term t1)
{
	bool flag = false;//no such term by efault
	trigPowers cur;
	cur.value = t1.value;
	cur.sinPower = 0;
	cur.cosPower = 0;

	if (t1.len != -1)// not C-case
	{
		cur.value /= t1.len;
		for (unsigned int i = 0; i < t1.len; i++)
			if (t1.ops[i] == 'z'){
				cur.cosPower++;
			}
			else {
				cur.sinPower++;
				cur.value *= Cos::getSign(t1.nums[i]);
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

void groundEnergy::addTermRotationZY2sublattices(int order, term t1)
{
	bool flag = false;//no such term by efault
	trigPowers cur;
	cur.value = t1.value;
	cur.sinPower = 0;
	cur.cosPower = 0;
	cur.complex_power = 0;

	if (t1.len != -1)// not C-case
	{
		cur.value /= t1.len;
		for (unsigned int i = 0; i < t1.len; i++)
			if (t1.ops[i] == 'z'){
				cur.cosPower++;
			}
			else if (t1.ops[i] == 'p') {
				cur.sinPower++;
				cur.value *= Cos::getSign(t1.nums[i]);
				cur.complex_power++;
				if (cur.complex_power == 2){
					cur.value *= -1;
					cur.complex_power = 0;
				}

			}
			else if (t1.ops[i] == 'm') {
				cur.sinPower++;
				cur.value *= -Cos::getSign(t1.nums[i]);
				cur.complex_power++;
				if (cur.complex_power == 2){
					cur.value *= -1;
					cur.complex_power = 0;
				}
			}
			else {
				std::cout << "WRONG ENERGY";
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

void groundEnergy::printTermRotation(std::ostream& out, int order)
{
	for (unsigned int i = 0; i < trigEnergy[order].size(); i++)
	{
		if (trigEnergy[order][i].complex_power != 0 && abs(trigEnergy[order][i].value) > 0.00000001)
			std::cout << "AAAAAAAAAAAAAAAA\n AAAAAAAAAAAAAAAAAAAAAAA\n AAAAAAAAAAAAAAAAAA\n !!!!!!!!!!!!!!!!";
		out << trigEnergy[order][i].value;
		switch (trigEnergy[order][i].cosPower){
			case 0: break;
			case 1: out << "*Cos[beta]"; break;
			default:out << "*Cos[beta]^" << trigEnergy[order][i].cosPower; break;
		} 
		switch (trigEnergy[order][i].sinPower){
			case 0: break;
			case 1: out << "*Sin[beta]"; break;
			default: out << "*Sin[beta]^" << trigEnergy[order][i].sinPower; break;
		}
		if (!ifNumerical){
			switch (trigEnergy[order][i].cosPower + trigEnergy[order][i].sinPower){
				case 0: break;
				case 1: out << "*Sz"; break;
				default: out << "*Sz^" << (trigEnergy[order][i].cosPower + trigEnergy[order][i].sinPower);
			}
		}
		else {
			switch (trigEnergy[order][i].cosPower + trigEnergy[order][i].sinPower){
				case 0: break;
				case 1: out << "*(" << spin << ")"; break;
				default: out << "*(" << spin << ")^" << (trigEnergy[order][i].cosPower + trigEnergy[order][i].sinPower);
			}
			
		}
			
		if (i != trigEnergy[order].size() - 1)
			out << "+";
	}
}
