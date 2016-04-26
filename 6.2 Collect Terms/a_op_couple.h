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
	size_t hash;
	bool is_hash_set;

	a_op_couple();

	bool operator==(const a_op_couple sec) const;

	bool operator<(const a_op_couple sec) const;

	void shift_in_elementary_cell()
	{
		switch (names[0])
		{
		case 'm': names[0] = 'M'; break;
		case 'p': names[0] = 'P'; break;
		case 'M': names[0] = 'm'; break;
		case 'P': names[0] = 'p'; break;
		}
		switch (names[1])
		{
		case 'm': names[1] = 'M'; break;
		case 'p': names[1] = 'P'; break;
		case 'M': names[1] = 'm'; break;
		case 'P': names[1] = 'p'; break;
		}
	}

	void printAterm(std::ofstream &F, int **m, int size, bool if_print_coeff = true)const;

	void printAtermTransfer(std::ofstream &F, int **m, int size, bool if_print_coeff = true)const;

	void printDoubleAterm(std::ofstream &F, int **m, int size, bool if_print_coeff = true)const;

	void check(); //dx set to be always greater or equal to zero

	void setHash(bool force = false);

	size_t getHash() const;

	void clear();

};
#endif //end of __A_OP_COUPLE_H_INCLUDED__