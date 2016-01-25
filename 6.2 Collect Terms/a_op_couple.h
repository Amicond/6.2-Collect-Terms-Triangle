#include "stdafx.h"
#ifndef __A_OP_COUPLE_H_INCLUDED__  
#define __A_OP_COUPLE_H_INCLUDED__

class a_op_couple{ //special class for bilinear terms
public:
	char names[2]; //plus or minus
	int dx,dy; //shift
	double coeff;
	int order;
	int i_power;

	a_op_couple();

	bool operator==(const a_op_couple sec) const;

	bool operator<(const a_op_couple sec) const;

	void printAterm(std::ofstream &F, int **m, int size, bool if_print_coeff = true)const;

	void check(); //dx set to be always greater or equal to zero

	void clear();

};
#endif //end of __A_OP_COUPLE_H_INCLUDED__