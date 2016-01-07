#include "stdafx.h"
#include "a_op.h"
#include "cos.h"

bool a_op::operator==(const a_op sec) const
{
	if (n != sec.n)
		return false;
	for (int i = 0; i < n; i++)
	{
		if (names[i] != sec.names[i])
			return false;
		if (node[i] != sec.node[i])
			return false;
	}
	return true;
}

bool a_op::operator<(const a_op sec) const
{
	if (n < sec.n) return true;
	else if (n > sec.n) return false;
	for (int i = 0; i < n; i++)	{
		if (names[i] < sec.names[i]) return true;
		else if (names[i] > sec.names[i]) return false;
	}
	for (int i = 0; i < n; i++)	{
		if (node[i] < sec.node[i]) return true;
		else if (node[i] > sec.node[i]) return false;
	}
	return false;
}

void a_op::printAterm(std::ofstream &F, int **m, int size, bool if_print_coeff)const
{
	if (n == 2)
	{
		if (if_print_coeff)
			F << coeff << "*";
		if (names[0] != names[1])
			F << "G";
		else if (names[0] == 'm')
			F << "F";
		else if (names[0] == 'p')
			F << "Fp";
		else
			F << "\n\n Strange ops \n\n";
		if (node[0] != node[1])
		{
			//if (node[0]!=0)
			//	F << "\n\n Strange nodes \n\n";
			int da1, db1, da2, db2;
			Cos::findCos( node[0], da1, db1);
			Cos::findCos( node[1], da2, db2);
			F << "*Cos[";
			switch ((da2 - da1)){
			case 0: break;
			case -1: F << "-ka"; break;
			case 1: F << "ka"; break;
			default: F << "ka*" << (da2 - da1);
			}
			switch ((db2 - db1)){
			case 0: break;
			case -1: F << "-kb"; break;
			case 1: F << "+kb"; break;
			default: F << "+kb*" << (db2 - db1);
			}
			F << "]";
		}
	}
	else
	{
		F << "\n\n Strange length \n\n";
	}
}

void a_op::clear()
{
	int n = 0;
	names.clear();
	node.clear();
	double coeff = 0;
	int order = 0;
}