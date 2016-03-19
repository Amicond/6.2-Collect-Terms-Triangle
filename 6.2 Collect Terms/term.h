#include "stdafx.h"
#ifndef __TERM_H_INCLUDE__
#define __TERM_H_INCLUDE__
class term
{
public:
	bool type; //if z-only or C then true(for energy), else false -correction
	int len;
	int order;
	char ops[10];
	int nums[10];
	double value;
	size_t hash;
	bool is_hash_set;
public:
	void decompose(std::string s, double val);
	
	void setOrder(int Ord);

	bool operator==(const term next) const;
	bool operator<(const term t2)const;
	
	void moveToZero();
	
	size_t getHash() const;
	void term::setHash();

};
#endif //__TERM_H_INCLUDE__