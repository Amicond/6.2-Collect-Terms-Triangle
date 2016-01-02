#include "stdafx.h"
#include "correction.h"

bool correction::operator == (correction c2)
{
	for (unsigned int i = 0; i < in.size(); i++)
	{
		if (in[i] != c2.in[i])
			return false;
	}
	if (out[0] != c2.out[0] || out[1] != c2.out[1])
		return false;
	return true;
}
void correction::add(Cos new_cos)
{
	std::vector<Cos>::iterator it;
	it = std::find(cs.begin(), cs.end(), new_cos);
	if (it != cs.end())
	{
		it->factor += new_cos.factor;
	}
	else
		cs.push_back(new_cos);
}
void correction::print(std::ofstream & outF, int num)
{
	if (num == 1)//for sz terms
	{
		outF << "Subscript[ap,k]*Subscript[a,k]*" << cs[0].factor;
	}

	if (num == 2)  //for pm terms
	{
		if (out[0] == 'p')
			outF << "Subscript[ap,k]*";
		else
			outF << "Subscript[a,k]*";

		if (out[0] == 'p'&&out[1] == 'p')
			outF << "Subscript[ap,-k]";
		if (out[0] == 'p'&&out[1] == 'm')
			outF << "Subscript[a,k]";
		if (out[0] == 'm'&&out[1] == 'p')
			outF << "Subscript[ap,k]";
		if (out[0] == 'm'&&out[1] == 'm')
			outF << "Subscript[a,-k]";

		outF << "*(";
		for (unsigned int i = 0; i < cs.size(); i++)
		{
			outF << cs[i].factor << "*Cos[";
			if (out[0] == 'p')
				outF << "Subscript[k,a]*" << cs[i].ka[0] << "+Subscript[k,b]*" << cs[i].kb[0] << "+";
			else
				outF << "-Subscript[k,a]*" << cs[i].ka[0] << "-Subscript[k,b]*" << cs[i].kb[0] << "+";
			///////////////////////////////////////////////////////////////////////////////////////////////
			// cos
			///////////////////////////////////////////////////////////////////////////////////////////////

			if (out[0] == 'p')
				outF << "-Subscript[k,a]*" << cs[i].ka[1] << "-Subscript[k,b]*" << cs[i].kb[1];
			else
				outF << "Subscript[k,a]*" << cs[i].ka[1] << "+Subscript[k,b]*" << cs[i].kb[1];

			outF << "]";
			if (i != cs.size() - 1)
			{
				outF << "+";
			}
		}
		outF << ")";

	}
	if (num == 3)  //for triple terms
	{
		if (in[0] == 'p')
			outF << "Subscript[ap,k1]*";
		else
			outF << "Subscript[a,k1]*";
		if (out[0] == 'p')
			outF << "Subscript[ap,k2]*";
		else
			outF << "Subscript[a,k2]*";
		if (out[1] == 'p')
			outF << "Subscript[ap,k3]*";
		else
			outF << "Subscript[a,k3]*";
		//temp
		//char c=
		//end temp
		outF << "DiracDelta[" << ((in[0] == 'p') ? "k1" : "-k1") << "+" << (out[0] == 'p' ? "k2" : "-k2") << "+" << (out[1] == 'p' ? "k3" : "-k3") << "]*(";
		for (unsigned int i = 0; i < cs.size(); i++)
		{
			outF << cs[i].factor << "*Cos[";
			if (in[0] == 'p')
				outF << "Subscript[k1,a]*" << cs[i].ka[0] << "+Subscript[k1,b]*" << cs[i].kb[0] << "+";
			else
				outF << "Subscript[k1,a]*" << -cs[i].ka[0] << "+Subscript[k1,b]*" << -cs[i].kb[0] << "+";
			if (out[0] == 'p')
				outF << "Subscript[k2,a]*" << cs[i].ka[1] << "+Subscript[k2,b]*" << cs[i].kb[1] << "+";
			else
				outF << "Subscript[k2,a]*" << -cs[i].ka[1] << "+Subscript[k2,b]*" << -cs[i].kb[1] << "+";
			if (out[1] == 'p')
				outF << "Subscript[k3,a]*" << cs[i].ka[1] << "+Subscript[k3,b]*" << cs[i].kb[1];
			else
				outF << "Subscript[k3,a]*" << -cs[i].ka[1] << "+Subscript[k3,b]*" << -cs[i].kb[1];
			outF << "]";
			if (i != cs.size() - 1)
			{
				outF << "+";
			}
		}
		outF << ")";

	}
	if (num == 4)  //for quatro terms
	{
		if (in[0] == 'p'&&in[1] == 'p')
			outF << "Subscript[Fp,k]";
		if (in[0] == 'm'&&in[1] == 'm')
			outF << "Subscript[Fm,k]";
		if (in[0] == 'p'&&in[1] == 'm')
			outF << "Subscript[G,k]";
		if (in[0] == 'm'&&in[1] == 'p')
			outF << "Subscript[G,k]";
		outF << "*";
		if (out[0] == 'p'&&out[1] == 'p')
			outF << "Subscript[Fp,q]";
		if (out[0] == 'm'&&out[1] == 'm')
			outF << "Subscript[Fm,q]";
		if ((out[0] == 'p'&&out[1] == 'm') || (out[0] == 'm'&&out[1] == 'p'))
			outF << "Subscript[G,q]";
		outF << "*(";
		for (unsigned int i = 0; i < cs.size(); i++)
		{
			outF << cs[i].factor << "*Cos[";
			outF << "ka*" << cs[i].ka[0] << +"+kb*" << cs[i].kb[0] << "+";
			if (in[1] == in[0])
				outF << "ka*" << cs[i].ka[1] << +"+kb*" << cs[i].kb[1];
			else
				outF << "ka*" << -cs[i].ka[1] << +"+kb*" << -cs[i].kb[1];
			outF << "+";
			outF << "qa*" << cs[i].ka[2] << +"+qb*" << cs[i].kb[2];
			outF << "+";
			if (out[1] == out[2])
				outF << "qa*" << cs[i].ka[3] << +"+qb*" << cs[i].kb[3];
			else
				outF << "qa*" << -cs[i].ka[3] << +"+qb*" << -cs[i].kb[3];
			outF << "]";
			if (i != cs.size() - 1)
			{
				outF << "+";
			}
		}
		outF << ")";
	}


}