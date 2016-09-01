#include "stdafx.h"
#ifndef __TERM_STORAGE__
#define __A_OP_XZ_ROTATE_STORAGE__
#include "term.h"
#include "cos.h"

namespace std{
	template <>
	struct hash < term >
	{
		std::size_t operator()(const term & k) const; 
	};
}

class term_storage{

	std::unordered_map<term,double> terms_storage;
	
	//Methods
	std::vector<std::pair<int, int>> generate_pairs(int n);
public:

	void clear_terms();

	void add_term(const term &t1);

	void print_terms(std::ofstream &F);
};
#endif
