#include "stdafx.h"
#ifndef __COS_H__
#define __COS_H__
class Cos
{
public:
	static int **m; //matrix with node-numbers
	static int size; //size of the matrix
	double factor;
	std::vector<int> ka, kb;
	static std::vector<std::pair<int, int>> coords;
	bool operator ==(Cos c2);
	static void set(int **M, int  Size);
	static void findCos( int n, int &da, int &db);
	static void findArbitraryCos(int n1, int n2,int &da, int &db);
	static int getSign(int num);
	static int getSignPiZero(int num); //different signs on vertical rows
	static int getSignZeroPi(int num); //different signs on horizontal rows
	static void printArbitraryCos(std::ofstream &F, int n1, int n2);
	static void newNodeAfterShift(int inp, int da, int db, int &out);
};
#endif //__COS_H__