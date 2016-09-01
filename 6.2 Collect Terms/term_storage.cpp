#include "stdafx.h"
#include "term.h"
#include "term_storage.h"


std::size_t std::hash<term>::operator()(const term & k) const {
	if (k.is_hash_set)
		return k.hash;
	else
		return k.getHash();
}

void term_storage::clear_terms()
{
	terms_storage.clear();
}

void term_storage::add_term(const term &t1)
{
	if (t1.len < 2) return;

	std::unordered_map<term, double>::iterator it;
	std::vector<std::pair<int, int>> pairs = generate_pairs(t1.len);
	for (int i = 0; i < t1.len; i++)
		if (t1.ops[i] != 'z') 
			return;
	for (auto &elem : pairs)
	{
		term cur_t;
		cur_t.len = 2;
		cur_t.value = t1.value;
		cur_t.order = t1.order;
		cur_t.ops[0] = t1.ops[elem.first];
		cur_t.ops[1] = t1.ops[elem.second];
		cur_t.nums[0] = t1.nums[elem.first];
		cur_t.nums[1] = t1.nums[elem.second];
		cur_t.moveToZero();
		cur_t.setHash();
		for (int i = 0; i < t1.len - 2; i++)
			cur_t.value *= -0.5;
		it = terms_storage.find(cur_t);
		if (it != terms_storage.end())
			it->second += cur_t.value;
		else
			terms_storage.insert({ cur_t, cur_t.value });
	}
}


std::vector<std::pair<int, int>> term_storage::generate_pairs(int n)
{
	//√енерирует все возможные пары дл€ данного терма
	std::vector<std::pair<int, int>> pairs;
	for (int i = 0; i < n; i++)	{
		for (int j = i + 1; j < n; j++) {
			pairs.push_back(std::pair<int, int>::pair(i, j));
		}
	}
	return pairs;
}

void term_storage::print_terms(std::ofstream &F) {
	if (terms_storage.size() == 0)
		F << "0";
	else {
		std::unordered_map<term, double>::iterator it;
		it = terms_storage.begin();
		while (it != terms_storage.end())
		{
			F << it->second << "*";
			Cos::printArbitraryCos(F, it->first.nums[0], it->first.nums[1]);
			it++;
			if (it != terms_storage.end())
				F << "+";
		}
	}
}