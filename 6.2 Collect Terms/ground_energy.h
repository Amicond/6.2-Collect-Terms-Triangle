#include "stdafx.h"
#ifndef __GROUND_ENERGY_H__
#define __GROUND_ENERGY_H__
#include "term.h"
#include "pseudospin_types.h"


const int N = 8; //max order


//class for several types of pseudospins
class groundEnergy
{
	struct trigPowers //for coeeficients in case of rotation
	{
		int cosPower[Types];
		int sinPower[Types];

		double value;

		bool operator==(const trigPowers& right)const;
		
		trigPowers();

		void add_trig(int type, int power, int plaquet_type);
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
	

	//void addTermRotationZeroPi2Angles(int order, term t1);

	//void addTermRotationPiZero2Angles(int order, term t1);

	//void addTermRotation2sublattices(int order, term t1);

	void addTriangleTermRotation(int order, term t1, int zeroType);
	
	void printTermRotation(std::ostream& out, int order, bool two_angels=false);
	
};

#endif //__GROUND_ENERGY_H__