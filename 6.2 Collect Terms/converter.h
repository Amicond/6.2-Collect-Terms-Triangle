#include "stdafx.h"
#ifndef __CONVERTER_H__ 
#define __CONVERTER_H__
#include "correction.h"
#include "a_op.h"

//forward declaration
class term;

//arrays
extern int mask4[6][4];
extern int mask3[3][3];
extern int out_pair[6][3];

extern std::string green_func[4];
extern std::string momenta_names[3];
//class defenition
class converter
{
	int **m;
	int matrix_size;
	std::vector<term> shorterTerms;
	std::vector<a_op> a_ops_l;
	//vector<a_op_rot> a_ops_l_rot;
	std::vector<correction> cors;
	double factor;
	int a_amount;

private:
	void insertShortTerm(term cur);
	void decomposeTerm1(term in);
	void decomposeTerm2(term in);
	void decomposeTerm3(term in);
	void decomposeTerm4(term in);
public:
	void set(double F, int A_amount, int **M, int Matrix_size);
	void decomposeTerm(term in);
	void convertToAop();
	void convertToCorrections();
	void set_signs(char op1, char op2, int &s1, int &s2, int &type);
	void print_3_momenta(std::ostringstream& outF, int signs1[], char type, int k);
	void combine(std::ofstream &outF);
	void clearTerms();
	bool PrintAll(std::ofstream &F);
	bool PrintAllAop(std::ofstream &F);
	
};

#endif //end if __CONVERTER_H__