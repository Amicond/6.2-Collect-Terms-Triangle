#include "stdafx.h"
#ifndef __GROUND_ENERGY_H__
#define __GROUND_ENERGY_H__
#include "term.h"


const int N = 8; //max order


//class for several types of pseudospins
class groundEnergy
{
	static const int Types = 3;
	struct trigPowers //for coeeficients in case of rotation
	{
		int cosPower[;
		int sinPower;
		int cosPower2;
		int sinPower2;
		double value;
		int complex_power;//for complex i, in case of ZY rotation

		bool operator==(const trigPowers& right)const;
		
		trigPowers();

		void add_trig(int type, int power, bool row_type);
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

	void addTermRotationZeroPi2Angles(int order, term t1);

	void addTermRotationPiZero2Angles(int order, term t1);

	void addTermRotation2sublattices(int order, term t1);

	void addTermRotationZY2sublattices(int order, term t1);
	
	void printTermRotation(std::ostream& out, int order, bool two_angels=false);
	
};

#endif //__GROUND_ENERGY_H__