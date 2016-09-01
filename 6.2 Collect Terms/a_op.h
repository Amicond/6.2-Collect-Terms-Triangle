#include "stdafx.h" 
#ifndef __A_OP_H_INCLUDED__  
#define __A_OP_H_INCLUDED__
class a_op
{
public:
	int n;//general length
	std::vector<char> names; //plus or minus
	std::vector<int> node; //to calc shift
	double coeff;
	int order;
	bool operator==(const a_op sec) const;

	bool operator<(const a_op sec) const;
	
	void printAterm(std::ofstream &F, int **m, int size, bool if_print_coeff = true)const;
	

	void clear();

};
#endif //end of __A_OP_H_INCLUDED__