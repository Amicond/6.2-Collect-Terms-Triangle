#include "stdafx.h"
#ifndef __COS_H__
#define __COS_H__
class Cos
{
public:
	double factor;
	std::vector<int> ka, kb;
	bool operator ==(Cos c2);
	static void Cos::findCos(int **m, int  size, int n, int &da, int &db);
};
#endif //__COS_H__