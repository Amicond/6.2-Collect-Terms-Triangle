#include "stdafx.h"
#include "a_op_couple.h"
#include "cos.h"



a_op_couple::a_op_couple()
{
	i_power = 0;
	is_hash_set = false;
}

bool a_op_couple::operator==(const a_op_couple sec) const
{
	if (names[0] != sec.names[0])
		return false;
	if (names[1] != sec.names[1])
		return false;
	if (dx != sec.dx)
		return false;
	if (dy != sec.dy)
		return false;
	if (i_power != sec.i_power)
		return false;
	return true;
}

size_t a_op_couple::getHash() const
{
	if (!is_hash_set){
		std::ostringstream out;
		out << names[0] << " " << names[1] << " ";
		out << dx << " " << dy << " ";
		out << i_power;
		return (std::hash<std::string>()(out.str()));
	}
	else
		return hash;
}

void a_op_couple::setHash()
{
	if (!is_hash_set){
		is_hash_set = true;
		hash = getHash();
	}
}

void a_op_couple::check()
{
	if (dx < 0)//always positive dx
	{
		dx = -dx;
		dy = -dy;
	}
	else if (dx==0&&dy<0)
	{
		dy = -dy;
	}
	if (names[0] != names[1] && names[0] == 'm') // G-function case for excitations
	{
		names[0] = 'p';
		names[1] = 'm';
	}
}

bool a_op_couple::operator<(const a_op_couple sec) const
{
	
	for (int i = 0; i < 2; i++)	{
		if (names[i] < sec.names[i]) return true;
		else if (names[i] > sec.names[i]) return false;
	}
	
	if (dx < sec.dx) return true;
	else if (dx > sec.dx) return false;

	if (dy < sec.dy) return true;
	else  return false;
}

void a_op_couple::printAterm(std::ofstream &F, int **m, int size, bool if_print_coeff)const
{
	
	if (if_print_coeff)
		F << coeff << "*";
	if (i_power != 0)
		F << "Sqrt[-1]^(" << i_power << ")*";
	if (names[0] != names[1])
		F << "G";
	else if (names[0] == 'm')
		F << "F";
	else if (names[0] == 'p')
		F << "Fp";
	else
		F << "\n\n Strange ops \n\n";
	if (dx != 0 || dy != 0)
	{
		//if (node[0]!=0)
		//	F << "\n\n Strange nodes \n\n";
		F << "*Cos[";
		switch (dx){
		case 0: break;
		case 1: F << "ka"; break;
		default: F << "ka*" << dx;
		}
		switch (dy){
		case 0: break;
		case -1: F << "-kb"; break;
		case 1: F << "+kb"; break;
		default: F << "+kb*" << (dy);
		}
		F << "]";
	}
	
}

void a_op_couple::printAtermTransfer(std::ofstream &F, int **m, int size, bool if_print_coeff) const
{
	if (if_print_coeff)
		F << coeff << "*";
	// test on complexity for not typical rotations
	if (i_power != 0)
		F << "Sqrt[-1]^(" << i_power << ")*";

	// names output (G=(a+_k)(a_k) F=(a_k)(a_-k) Fp=(a+_k)(a+_-k) GG=(a+_(k+k0))(a_k) FF=(a_k)(a_(-k+k0)) FFp=(a+_k)(a+_(-k+k0))
	//
	if (names[0] =='p'&& names[1]=='m')
		F << "G";
	else if (names[0] == 'm'&&names[1]=='m')
		F << "F";
	else if (names[0] == 'p'&&names[1] == 'p')
		F << "Fp";
	else if (names[0] == 'P'&&names[1] == 'm')
		F << "GG";
	else if (names[0] == 'P'&&names[1] == 'p')
		F << "FFp";
	else if (names[0] == 'M'&&names[1] == 'm')
		F << "FF";
	else 
		F << "\n\n Strange ops \n\n";

	//print shift cosine
	if (dx != 0 || dy != 0)
	{
		//if (node[0]!=0)
		//	F << "\n\n Strange nodes \n\n";
		F << "*Cos[";
		switch (dx) {
			case 0: break;
			case 1: F << "ka"; break;
			default: F << "ka*" << dx;
		}
		switch (dy) {
			case 0: break;
			case -1: F << "-kb"; break;
			case 1: F << "+kb"; break;
			default: F << "+kb*" << (dy);
		}
		F << "]";
	}

}

void a_op_couple::clear()
{
	coeff = 0;
	order = 0;
	i_power = 0;
}