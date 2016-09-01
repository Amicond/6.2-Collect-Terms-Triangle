#include "stdafx.h"
#include "pseudospin_types.h"

const std::string pseudospin_types::angleNames[Types] = { "alpha","beta","gamma" };

int ** pseudospin_types::matrix;
int pseudospin_types::NNN;
std::string pseudospin_types::getAngleName(int i)
{
	return angleNames[i];
}

int pseudospin_types::getType(int n, int zeroType)//проверить
{
	int x, y;
	int d = size();
	for (int i = 0; i < d;i++)
		for (int j = 0; j < d; j++)
			if(matrix[i][j]==n)
			{ 
				x = j;
				y = i;
				break;
			}
	x = x - d / 2;
	y = y - d / 2;
	return ((zeroType + (y - x)) % 3 + 3) % 3;;
}

int pseudospin_types::size()
{
	return std::max(1 + ((NNN) / 2) * 2, 1 + (NNN - 2) * 2);//
}

void pseudospin_types::setOrder(int NNN)
{
	pseudospin_types::NNN = NNN;
}


void pseudospin_types::fillMatrix()
{
	int dir = 0;
	int d = size();

	matrix = new int*[d];
	for (int i = 0; i<d; i++)
		matrix[i] = new int[d];


	int s = d / 2;
	int cx, cy, ct, cb, cl, cr;
	dir = 0;//0- направо, 1-вверх, 2- налево, 3- вниз
	ct = cb = cl = cr = cx = cy = s;
	for (int curn = 0; curn<d*d; curn++)
	{
		matrix[cy][cx] = curn;
		switch (dir)
		{
		case 0:
			if (cx == cr)
			{
				cr++;
				cl--;
				cx++;
				dir = 1;
			}
			else
			{
				cx++;
			}
			break;
		case 1:
			if (cy == ct)
			{
				ct--;
				cb++;
				cy--;
				dir = 2;
			}
			else
			{
				cy--;
			}
			break;
		case 2:
		{
			cx--;
			if (cx == cl)
			{
				dir = 3;
			}
		}
		break;
		case 3:
		{
			cy++;
			if (cy == cb)
			{
				dir = 0;
			}
		}
		break;
		}
	}
}
