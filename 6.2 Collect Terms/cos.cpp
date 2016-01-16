#include "stdafx.h"
#include "cos.h"

int ** Cos::m = NULL;
int Cos::size = -1;
std::vector<std::pair<int, int>> Cos::coords;

void Cos::findCos( int n, int &da, int &db)
{
	/*for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			if (m[i][j] == n)
			{
				da = j;
				db = i;
				break;
			}
		}*/
	da = coords[n].first;
	db = coords[n].second;
	da = da - size / 2;
	db = -db + size / 2;
}

void Cos::findArbitraryCos(int n1, int n2, int &da, int &db)
{
	/*for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) {
			if (m[i][j] == n1) {
				da = j;
				db = i;
				break;
			}
		}
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) {
			if (m[i][j] == n2) {
				da -= j;
				db -= i;
				break;
			}
		}*/
	da = coords[n1].first - coords[n2].first;
	db = coords[n1].second - coords[n2].second;
	//if (da1 != da || db1 != db) std::cout << "\n" << n1 << " " << n2 << "\n";
	if (da < 0){
		da = -da;
		db = -db;
	}
}

void Cos::set(int **M, int  Size){
	Cos::m = M;
	Cos::size = Size;
	for (int i = 0; i < Cos::size*Cos::size; i++) {
		coords.push_back(std::pair<int,int>({ 0, 0 }));
	}
	for (int i = 0; i < Cos::size; i++) {
		for (int j = 0; j < Cos::size; j++) {
			coords[m[i][j]].first = j;
			coords[m[i][j]].second = i;
		}
	}
	//for tets
	Cos::size = Size;
	//end test
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

int Cos::getSign(int num)
{
	int da, db;
	findCos(num, da, db);
	if ((da + db) % 2 == 0) return 1;
	else return -1;
}