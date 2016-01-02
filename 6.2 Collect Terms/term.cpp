#include "stdafx.h"

#include "term.h"
void term::decompose(std::string s, double val)
{
	int temp = s.find_first_not_of("SJpmz");
	//temp
	//int nn = s.find_first_of("pm");
	//end temp
	type = s.find_first_of("pm") == std::string::npos ? true : false;
	len = temp - 1;
	for (int i = 0; i < len; i++)
		ops[i] = s[i + 1];
	std::istringstream is;
	char tc;
	is.str(s.substr(temp));
	for (int i = 0; i <len; i++)
	{
		is >> nums[i];
		is >> tc;
	}
	value = val;

}
void term::setOrder(int Ord)
{
	order = Ord;
}

bool term::operator==(const term next) const
{
	if (len != next.len)
		return false;
	for (int i = 0; i < len; i++)
	{
		if (ops[i] != next.ops[i])
			return false;
		if (nums[i] != next.nums[i])
			return false;
	}
	return true;
}
bool term::operator<(const term t2)const
{
	if (order < t2.order) return true;
	else if (order > t2.order) return false;
	else{
		if (len < t2.len) return true;
		else if (len > t2.len) return false;
		else{
			for ( int i = 0; i < len; i++){
				if (ops[i] < t2.ops[i]) return true;
				else if (ops[i] > t2.ops[i])return false;
			}
			for ( int i = 0; i < len; i++){
				if (nums[i] < t2.nums[i]) return true;
				else if (nums[i] > t2.nums[i]) return false;
			}
			return false;
		}
	}
}
