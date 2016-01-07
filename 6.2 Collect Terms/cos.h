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
	bool operator ==(Cos c2);
	static void set(int **M, int  Size);
	static void Cos::findCos( int n, int &da, int &db);
	static void Cos::findArbitraryCos(int n1, int n2,int &da, int &db);
	static int getSign(int num);
	
};
#endif //__COS_H__