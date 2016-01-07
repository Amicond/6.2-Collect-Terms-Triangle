#include "stdafx.h"
#include "converter.h"
#include "term.h"
#include "a_op.h"


//arrays
int mask4[6][4] = { { 1, 1, 0, 0 }, { 1, 0, 1, 0 }, { 1, 0, 0, 1 }, { 0, 1, 1, 0 }, { 0, 1, 0, 1 }, { 0, 0, 1, 1 } };
int mask3[3][3] = { { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 } };
int out_pair[6][3] = { { 0, 1, 2 }, { 1, 0, 2 }, { 1, 2, 0 }, { 0, 2, 1 }, { 2, 0, 1 }, { 2, 1, 0 } };

std::string green_func[4] = { "G", "GP", "F", "FP" };
std::string momenta_names[3] = { "k", "k1", "k2" };

//methods

void converter::set(double F, int A_amount, int **M, int Matrix_size)
{
	factor = F;
	a_amount = A_amount;
	m = M;
	matrix_size = Matrix_size;
}

void converter::decomposeTerm(term in)
{
	if (a_amount == 1)
		decomposeTerm1(in);
	if (a_amount == 2)
	{
		decomposeTerm2(in);
	}
	if (a_amount == 3)
		decomposeTerm3(in);
	if (a_amount == 4)
		decomposeTerm4(in);
}


void converter::insertShortTerm(term cur)
{
	std::vector<term>::iterator it;
	it = find(shorterTerms.begin(), shorterTerms.end(), cur);
	if (it != shorterTerms.end())
		it->value += cur.value;
	else
		shorterTerms.push_back(cur);
}
void converter::decomposeTerm1(term in)
{
	term cur;
	cur.len = 1;
	double coeff = 1;
	for (unsigned int i = 0; i < in.len - 1; i++)
		coeff *= factor;
	cur.order = in.order;
	cur.value = in.value*coeff;
	cur.ops[0] = in.ops[0];
	cur.nums[0] = 0;
	insertShortTerm(cur);

}
void converter::decomposeTerm2(term in)
{
	//two term types:
	//pm+z and z-only
	//


	int pm = 0;
	int pmInd[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	int index = 0;
	for (unsigned int i = 0; i < in.len; i++) //ищем первый pm оператор
		if (in.ops[i] != 'z')
		{
			//debug out
			//if (index > 1)
			//cout << "Alarm!!!!!!!!!!!!!!!!\n";
			//end debug out
			pm++;
			pmInd[index++] = i;
		}
	term cur;
	double coeff = 1;
	if (pm == 2)
	{
		for (unsigned int i = 0; i < in.len - 2; i++)
			coeff *= factor;
		cur.len = 2;
		cur.order = in.order;
		cur.value = in.value*coeff / in.len;
		cur.ops[0] = in.ops[pmInd[0]];
		cur.nums[0] = in.nums[pmInd[0]];

		cur.nums[1] = in.nums[pmInd[1]];
		cur.ops[1] = in.ops[pmInd[1]];

		insertShortTerm(cur);

	}
	else if (pm == 0)
	{
		if (in.len)//if not C-term
		{
			for (unsigned int i = 0; i < in.len - 1; i++)//change all except one to factor
				coeff *= factor;
			cur.len = 1;						//1 Sz operator
			cur.order = in.order;
			cur.value = in.value*coeff;

			cur.ops[0] = in.ops[0];
			cur.nums[0] = in.nums[0];

			insertShortTerm(cur);
		}

	}

}
//TODO проверить делитель
void converter::decomposeTerm3(term in)
{
	term cur;
	cur.len = -1;
	double coeff = 1;

	int pm = 0;
	for (unsigned int i = 0; i < in.len; i++)
		if (in.ops[i] != 'z')
			pm++;
	if (pm == 1)
	{
		for (unsigned int i = 0; i < in.len - 2; i++)
			coeff *= factor;
		cur.len = 2;
		cur.order = in.order;
		cur.value = in.value*coeff / in.len;
		cur.ops[0] = in.ops[0];
		cur.nums[0] = in.nums[0];
		if (in.ops[0] == 'z')
		{
			for (unsigned int i = 1; i < in.len; i++) //ищем первый pm оператор
				if (in.ops[i] != 'z')
				{
					cur.nums[1] = in.nums[i];
					cur.ops[1] = in.ops[i];
					insertShortTerm(cur);
					break;
				}
		}
		else //выбираем все доступные z операторы
		{
			for (unsigned int i = 1; i < in.len; i++)
			{
				cur.ops[1] = in.ops[1];
				cur.nums[1] = in.nums[1];
				insertShortTerm(cur);
			}
		}
	}
	if (pm == 3)
	{
		for (unsigned int i = 0; i < in.len - 3; i++)
			coeff *= factor;
		cur.len = 3;
		cur.order = in.order;
		cur.value = in.value*coeff / in.len;
		int index = 0;
		for (unsigned int i = 0; i < in.len; i++)
		{
			if (in.ops[i] != 'z')
			{
				cur.ops[index] = in.ops[i];
				cur.nums[index] = in.nums[i];
				index++;
			}
		}
		insertShortTerm(cur);
	}

}
void converter::decomposeTerm4(term in)
{
	term cur;
	cur.len = -1;
	double coeff = 1;

	int pm = 0;
	for (unsigned int i = 0; i < in.len; i++)
		if (in.ops[i] != 'z')
			pm++;
	if (pm == 0)
	{
		for (unsigned int i = 0; i < in.len - 2; i++)
			coeff *= factor;
		cur.len = 2;
		cur.value = coeff*in.value / in.len;
		cur.order = in.order;
		for (unsigned int i = 1; i < in.len; i++)//заменяем все подряд sz, поочереди оставляем только i
		{

			cur.ops[0] = 'z';
			cur.ops[1] = 'z';
			cur.nums[0] = in.nums[0];
			cur.nums[1] = in.nums[i];

			//вставляем
			insertShortTerm(cur);
		}
	}
	if (pm == 2)
	{
		for (unsigned int i = 0; i < in.len - 3; i++)
			coeff *= factor;
		cur.len = 3;
		cur.value = coeff*in.value / in.len;
		cur.order = in.order;
		cur.nums[0] = in.nums[0];
		cur.ops[0] = in.ops[0];
		if (cur.ops[0] == 'z')//оставляем только pm-члены из оставшихся
		{
			int index = 1;
			for (unsigned int i = 1; i < in.len; i++)
			{
				if (in.ops[i] != 'z')
				{
					cur.ops[index] = in.ops[i];
					cur.nums[index] = in.nums[i];
					index++;
				}
			}
			//вставляем
			insertShortTerm(cur);
		}
		else //по очереди выбираем все z-члены
		{
			unsigned int index = 0;
			for (unsigned int i = 1; i < in.len; i++) //ищем второй pm оператор
				if (in.ops[i] != 'z')
				{
					index = i;
					break;
				}
			for (unsigned int i = 1; i < in.len; i++)
			{
				if (index < i)
				{
					cur.ops[1] = in.ops[index];
					cur.ops[2] = in.ops[i];
					cur.nums[1] = in.nums[index];
					cur.nums[2] = in.nums[i];
				}
				if (index > i)
				{
					cur.ops[1] = in.ops[i];
					cur.ops[2] = in.ops[index];
					cur.nums[1] = in.nums[i];
					cur.nums[2] = in.nums[index];
				}
				//вставляем
				if (i != index)
					insertShortTerm(cur);
			}
		}

	}

	if (pm == 4)
	{
		for (unsigned int i = 0; i < in.len - 4; i++)
			coeff *= factor;
		cur.len = 4;
		cur.order = in.order;
		cur.value = in.value*coeff / in.len;
		int index = 0;
		for (unsigned int i = 0; i < in.len; i++)
		{
			if (in.ops[i] != 'z')
			{
				cur.ops[index] = in.ops[i];
				cur.nums[index] = in.nums[i];
				index++;
			}
		}
		//вставляем
		insertShortTerm(cur);
	}
}

void converter::convertToAop()
{
	int pm;
	a_op cur;
	for (unsigned int i = 0; i < shorterTerms.size(); i++)
	{
		pm = 0;
		cur.names.clear();
		cur.node.clear();
		cur.n = 0;
		cur.coeff = shorterTerms[i].value;
		cur.order = shorterTerms[i].order;
		for (unsigned int j = 0; j < shorterTerms[i].len; j++)
		{
			if (shorterTerms[i].ops[j] == 'z')
			{
				cur.names.push_back('p');
				cur.names.push_back('m');
				cur.node.push_back(shorterTerms[i].nums[j]);
				cur.node.push_back(shorterTerms[i].nums[j]);
				cur.n += 2;
				if (factor > 0)//если sz>0, то надо брать знак минус
					cur.coeff *= -1;
			}
			else//pm-case
			{
				cur.names.push_back(shorterTerms[i].ops[j]);
				cur.node.push_back(shorterTerms[i].nums[j]);
				cur.n++;
			}
		}

		//vector<a_op>::iterator it = find(a_ops_l.begin(), a_ops_l.end(), cur);
		//if (it != a_ops_l.end())

		//start debug
		if (cur.node.size() >2)
		{
			std::cout << "Too long!!!!!!!!! test\n";
		}
		//end debug*/


		a_ops_l.push_back(cur);
	}
}

void converter::convertToCorrections()
{
	correction cur;
	int index, index2;
	if (a_amount == 1)//для Sz
	{
		for (unsigned int i = 0; i < a_ops_l.size(); i++)
		{
			cur.cs.clear();
			cur.out[0] = a_ops_l[i].names[0];
			cur.out[1] = a_ops_l[i].names[1];
			Cos cur_cos;
			cur_cos.factor = a_ops_l[i].coeff;

			cur_cos.ka.push_back(0);
			cur_cos.kb.push_back(0);

			cur_cos.ka.push_back(0);
			cur_cos.kb.push_back(0);
			std::vector<correction>::iterator it = find(cors.begin(), cors.end(), cur);
			if (it != cors.end())
			{
				it->cs.push_back(cur_cos);
			}
			else
			{
				cur.cs.push_back(cur_cos);
				cors.push_back(cur);
			}
		}
	}
	if (a_amount == 2)
	{
		for (unsigned int i = 0; i < a_ops_l.size(); i++)
		{
			cur.cs.clear();
			Cos cur_cos;
			cur_cos.factor = a_ops_l[i].coeff;
			int da, db;
			cur.out[0] = a_ops_l[i].names[0];
			cur.out[1] = a_ops_l[i].names[1];

			Cos::findCos( a_ops_l[i].node[0], da, db);
			cur_cos.ka.push_back(da);
			cur_cos.kb.push_back(db);
			Cos::findCos( a_ops_l[i].node[1], da, db);
			cur_cos.ka.push_back(da);
			cur_cos.kb.push_back(db);


			//end eval
			std::vector<correction>::iterator it = find(cors.begin(), cors.end(), cur);
			if (it != cors.end())
			{
				it->cs.push_back(cur_cos);
			}
			else
			{
				cur.cs.push_back(cur_cos);
				cors.push_back(cur);
			}

		}
	}
	if (a_amount == 3)
	{
		for (unsigned int i = 0; i < a_ops_l.size(); i++)
		{

			for (int j = 0; j < 3; j++) //3 способа выбрать внешний опреатор
			{
				cur.cs.clear();
				cur.in.clear();
				index = 0;
				int nums[3];
				for (int k = 0; k < 3; k++)//назначаем какой оператор будет внешним
				{
					if (mask3[j][k] == 0)//вставляем в in
					{
						cur.in.push_back(a_ops_l[i].names[k]);
						nums[0] = k;
					}
					else
					{
						nums[1 + index] = k;
						cur.out[index++] = a_ops_l[i].names[k];
					}
				}
				//eval cos for cur term
				Cos cur_cos;
				cur_cos.factor = a_ops_l[i].coeff;
				int da, db;

				Cos::findCos( a_ops_l[i].node[nums[0]], da, db);
				cur_cos.ka.push_back(da);
				cur_cos.kb.push_back(db);
				Cos::findCos( a_ops_l[i].node[nums[1]], da, db);
				cur_cos.ka.push_back(da);
				cur_cos.kb.push_back(db);
				Cos::findCos( a_ops_l[i].node[nums[2]], da, db);
				cur_cos.ka.push_back(da);
				cur_cos.kb.push_back(db);

				//end eval
				std::vector<correction>::iterator it = find(cors.begin(), cors.end(), cur);
				if (it != cors.end())
				{
					it->cs.push_back(cur_cos);
				}
				else
				{
					cur.cs.push_back(cur_cos);
					cors.push_back(cur);
				}

			}
		}
	}
	if (a_amount == 4)
	{
		for (unsigned int i = 0; i < a_ops_l.size(); i++)
		{

			for (int j = 0; j < 6; j++)
			{
				index = 0;
				index2 = 0;
				cur.cs.clear();
				cur.in.clear();
				int nums[4];
				for (int k = 0; k < 4; k++)
				{
					if (mask4[j][k] == 0)
					{
						cur.in.push_back(a_ops_l[i].names[k]);
						nums[index2++] = k;
					}
					else
					{
						nums[2 + index] = k;
						cur.out[index++] = a_ops_l[i].names[k];
					}
				}
				//eval cos for cur term
				Cos cur_cos;
				cur_cos.factor = a_ops_l[i].coeff;
				int da, db;

				Cos::findCos( a_ops_l[i].node[nums[0]], da, db);
				cur_cos.ka.push_back(da);
				cur_cos.kb.push_back(db);
				Cos::findCos( a_ops_l[i].node[nums[1]], da, db);
				cur_cos.ka.push_back(da);
				cur_cos.kb.push_back(db);
				Cos::findCos( a_ops_l[i].node[nums[2]], da, db);
				cur_cos.ka.push_back(da);
				cur_cos.kb.push_back(db);
				Cos::findCos( a_ops_l[i].node[nums[3]], da, db);
				cur_cos.ka.push_back(da);
				cur_cos.kb.push_back(db);

				//end eval
				std::vector<correction>::iterator it = find(cors.begin(), cors.end(), cur);
				if (it != cors.end())
				{
					it->cs.push_back(cur_cos);
				}
				else
				{
					cur.cs.push_back(cur_cos);
					cors.push_back(cur);
				}
			}
		}

	}
}

void converter::set_signs(char op1, char op2, int &s1, int &s2, int &type)
{
	if (op1 == 'm'&&op2 == 'p') //G_k, оба знака плюс
	{
		s1 *= 1;
		s2 *= 1;
		type = 0;
	}
	if (op1 == 'p'&&op2 == 'm') //Gp_k, оба знака минус
	{
		s1 *= -1;
		s2 *= -1;
		type = 1;
	}
	if (op1 == 'm'&&op2 == 'm') //F_k, + -
	{
		s1 *= 1;
		s2 *= -1;
		type = 2;
	}
	if (op1 == 'p'&&op2 == 'p') //F_k, - +
	{
		s1 *= -1;
		s2 *= 1;
		type = 3;
	}

}

void converter::print_3_momenta(std::ostringstream& outF, int signs1[], char type, int k)
	//type - 1: type==' ' - no additoions, 2: type=='a' - a axis, 3: 'b'- axis
{
	if (signs1[(k + 2) % 3] == 1)
	{
		if (-signs1[k] == -1)
			outF << "-" << momenta_names[0];
		else
			outF << momenta_names[0];
		if (type == 'a' || type == 'b') outF << type;
		if (-signs1[(k + 1) % 3] == -1)
			outF << "-" << momenta_names[1];
		else
			outF << "+" << momenta_names[1];
		if (type == 'a' || type == 'b') outF << type;
	}
	else
	{
		if (signs1[k] == -1)
			outF << "-" << momenta_names[0];
		else
			outF << momenta_names[0];
		if (type == 'a' || type == 'b') outF << type;
		if (signs1[(k + 1) % 3] == -1)
			outF << "-" << momenta_names[1];
		else
			outF << "+" << momenta_names[1];
		if (type == 'a' || type == 'b') outF << type;
	}
}

void converter::combine(std::ofstream &outF)
{
	int signs1[3], signs2[3];
	int da, db;
	int type[3];
	bool if_empty, total_empty1, total_empty2;
	std::ostringstream cos2;
	std::ostringstream momenta_3;
	if (a_amount == 3)
	{
		for (unsigned int i = 0; i < a_ops_l.size(); i++)//перебираем всех кого ставим в начало
		{
			for (unsigned int j = 0; j<a_ops_l.size(); j++) //назначаем всех поочереди внешним
			{
				for (unsigned int k = 0; k < 3; k++) //перебираем входной импульс из первых трех операторов
				{
					for (unsigned int l = 0; l < 6; l++) //задаем маску соответсвия, элементы массива указывают с каким оператором из второй тройки нужно спаривать
					{
						for (int mm = 0; mm < 3; mm++) //задаем знаки перед каждым импульсом по типу операторов
						{
							if (a_ops_l[i].names[mm] == 'p')
								signs1[mm] = 1;
							else
								signs1[mm] = -1;
							if (a_ops_l[j].names[mm] == 'p')
								signs2[mm] = 1;
							else
								signs2[mm] = -1;
						}

						for (int mm = 0; mm < 3; mm++)
						{
							set_signs(a_ops_l[i].names[mm], a_ops_l[j].names[out_pair[l][mm]], signs1[mm], signs2[out_pair[l][mm]], type[mm]);
						}
						cos2.str("");
						outF << "+" << a_ops_l[i].coeff*a_ops_l[j].coeff;
						//выводим 2 первых функции Грина
						for (int mm = 0; mm < 2; mm++)
						{
							outF << "*" << green_func[type[(k + mm) % 3]] << "_" << momenta_names[mm];
						}
						//вычисляем третий импульс через первые два
						outF << "*" << green_func[type[(k + 2) % 3]] << "_(";
						momenta_3.str("");
						print_3_momenta(momenta_3, signs1, ' ', k);
						outF << momenta_3.str();
						outF << ")";

						//

						outF << "*DiracDelta[";
						for (int mm = 0; mm < 3; mm++)
						{
							outF << "+(" << signs1[(k + mm) % 3] << ")*" << momenta_names[mm] << "a";
						}
						outF << "]*DiracDelta[";
						for (int mm = 0; mm < 3; mm++)
						{
							outF << "+(" << signs1[(k + mm) % 3] << ")*" << momenta_names[mm] << "b";
						}
						/*outF << "]*DiracDelta[";
						for (int mm = 0; mm < 3; mm++)
						{
						outF << "+(" << signs1[(k + mm) % 3] << ")*" << momenta_names[mm];
						}*/

						outF << "]*Cos[";
						cos2 << "*Cos[";
						total_empty1 = true;
						total_empty2 = true;
						for (int mm = 0; mm<3; mm++)
						{
							Cos::findCos( a_ops_l[i].node[(k + mm) % 3], da, db);
							if (da != 0 || db != 0)
							{
								total_empty1 = false;
								if (mm != 0)outF << "+";
								if (signs1[(k + mm) % 3] == -1)//знак плюс
								{
									da *= -1;
									db *= -1;
								}
								if_empty = true;
								if (da != 0)
								{
									if (mm != 2)
										outF << momenta_names[mm] << "a*" << da;
									else
									{
										momenta_3.str("");
										print_3_momenta(momenta_3, signs1, 'a', k);
										outF << "(" << momenta_3.str() << ")*" << da;
									}
									if_empty = false;
								}
								if (db != 0)
								{
									if (!if_empty) outF << "+"; //есть член с da "+", нет члена с da, + не нужен

									outF << momenta_names[mm] << "b*" << db;
								}

							}


							Cos::findCos( a_ops_l[j].node[out_pair[l][(k + mm) % 3]], da, db);
							//cos2 << "+" << signs2[out_pair[l][(k + mm) % 3]] << "*(" << momenta_names[mm] << "a*" << da << "+" << momenta_names[mm] << "b*" << db << ")";
							if (da != 0 || db != 0)
							{
								total_empty2 = false;
								if (mm != 0) cos2 << "+";
								if (signs2[out_pair[l][(k + mm) % 3]] == -1)//знак плюс
								{
									da *= -1;
									db *= -1;
								}
								if_empty = true;
								if (da != 0)
								{
									cos2 << momenta_names[mm] << "a*" << da;
									if_empty = false;
								}
								if (db != 0)
								{
									if (!if_empty)//есть член с da "+", нет члена с da, + не нужен
										cos2 << "+";
									cos2 << momenta_names[mm] << "b*" << db;
								}

							}
						}
						if (total_empty1)
							outF << "0";
						outF << "]";
						outF << cos2.str();
						if (total_empty2)
							outF << "0";
						outF << "]";

					}
				}
			}
		}
	}
}

void converter::clearTerms()
{
	shorterTerms.clear();
	a_ops_l.clear();
	cors.clear();
}

bool converter::PrintAll(std::ofstream &F)
{

	for (unsigned int i = 0; i < cors.size(); i++)
	{
		cors[i].print(F, a_amount);
		if (i != cors.size() - 1)
			F << "+";
	}
	if (cors.size() == 0)
		return false;
	else
		return true;
}

bool converter::PrintAllAop(std::ofstream &F)
{
	for (unsigned int i = 0; i < a_ops_l.size(); i++)
	{
		a_ops_l[i].printAterm(F, m, matrix_size);
		if (i != a_ops_l.size() - 1)
			F << "+";
	}
	if (a_ops_l.size() == 0)
		return false;
	else
		return true;
}

