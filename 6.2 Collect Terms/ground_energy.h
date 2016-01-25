#include "stdafx.h"
#ifndef __GROUND_ENERGY_H__
#define __GROUND_ENERGY_H__
#include "term.h"


const int N = 10; //max order

class groundEnergy
{
	struct trigPowers //for coeeficients in case of rotation
	{
		int cosPower;
		int sinPower;
		double value;
		int complex_power;//for complex i, in case of ZY rotation

		bool operator==(const trigPowers& right)const;
		
		trigPowers();
	};
	std::vector<double> energy;
	std::vector<trigPowers> trigEnergy[N];
	double spin;
	bool ifNumerical;

public:

	void set(double factor,bool IfNumerical);
	
	void addTerm(unsigned int order, term t1);

	void clearTerms();
	
	double returnE(unsigned int order);
	
	//rotation case
	void addTermRotation(int order, term t1);

	void addTermRotationAntiferromagnet(int order, term t1);

	void addTermRotationPiZero(int order, term t1);

	void addTermRotation2sublattices(int order, term t1);

	void addTermRotationZY2sublattices(int order, term t1);
	
	void printTermRotation(std::ostream& out, int order);
	
};

#endif //__GROUND_ENERGY_H__