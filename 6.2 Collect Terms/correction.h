#include "stdafx.h"
#ifndef __CORRECTIOSN__
#define __CORRECTIONS__
#include "cos.h"

class correction
{
public:
	std::vector<char> in;
	char out[2];
	std::vector<Cos> cs;
	bool operator ==(correction c2);
	
	void add(Cos new_cos);

	void print(std::ofstream & outF, int num);
};
#endif //__CORRECTIONS__