#include "stdafx.h"
#include "cos.h"

void Cos::findCos(int **m, int  size, int n, int &da, int &db)
{
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			if (m[i][j] == n)
			{
				da = j;
				db = i;
				break;
			}
		}
	da = da - size / 2;
	db = -db + size / 2;
}


bool Cos::operator ==(Cos c2)
{
	int s1, s2;
	s1 = ka.size();
	s2 = c2.ka.size();
	if (s1 != s2)
		return false;
	s1 = kb.size();
	s2 = c2.kb.size();
	if (s1 != s2)
		return false;
	for (unsigned int i = 0; i < ka.size(); i++)
	{
		if (ka[i] != c2.ka[i] || kb[i] != c2.kb[i])
			return false;
	}
	return true;
}