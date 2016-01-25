#include "stdafx.h"
#include "a_op_couple.h"
#include "a_op_xz_rotate_storage.h"



std::size_t std::hash<AopXZRotate>::operator()(const AopXZRotate& k) const{
	std::ostringstream out;
	out << k.aop_c.names[0] << " " << k.aop_c.names[1] << " ";
	out << k.aop_c.dx << " " << k.aop_c.dy << " ";
	out << k.sz_power << " ";
	out << k.trc.getPowers();
	return (std::hash<string>()(out.str()));
}

std::size_t std::hash<a_op_couple>::operator()(const a_op_couple& k) const {
	std::ostringstream out;
	out << k.names[0] << " " << k.names[1] << " ";
	out << k.dx << " " << k.dy << " ";
	out << k.i_power;
	return (std::hash<string>()(out.str()));
}


int AopXZRotateStorage::arr[4][2] = { { 0, 0 }, { 0, 1 }, { 1, 0 }, { 1, 1 } };

void AopXZRotateStorage::set(int **Matrix, int Matrix_size, bool Mode, double Angle_start, double Angle_finish, double Angle_step, double Angle_start_2, double Angle_finish_2, double Angle_step_2, double Sz)
{
	matrix = Matrix;
	matrix_size = Matrix_size;
	analytical_mode = Mode;

	if (analytical_mode == false) {

		angle_start = Angle_start;
		angle_step = Angle_step;
		angle_finish = Angle_finish;

		angle_start_2 = Angle_start_2;
		angle_step_2 = Angle_step_2;
		angle_finish_2 = Angle_finish_2;
		
		double angle_cur = angle_start;
		while (angle_cur <= angle_finish) {
			storage_numerical.push_back(std::unordered_map<a_op_couple, double>());
			angle_cur += Angle_step;
		}
	}
	sz = Sz;
}

void AopXZRotateStorage::clearTerms()
{
	storage.clear();
	for (unsigned int i = 0; i < storage_numerical.size(); i++) {
		storage_numerical[i].clear();
	}
	storage_numerical.clear();
}

std::vector<std::pair<int, int>> AopXZRotateStorage::generate_pairs(int n) {
	//√енерирует все возможные пары дл€ данного терма
	std::vector<std::pair<int, int>> pairs;
	for (int i = 0; i < n; i++)	{
		for (int j = i + 1; j < n; j++) {
			pairs.push_back(std::pair<int, int>::pair(i, j));
		}
	}
	return pairs;
}

void AopXZRotateStorage::add_numerical(AopXZRotate &current) {
	double angle_cur = angle_start;
	std::vector<std::unordered_map<a_op_couple, double>>::iterator it_vec = storage_numerical.begin();
	std::unordered_map<a_op_couple, double>::iterator it_map;
	double coeff_correction;
	while (angle_cur < angle_finish) {
		//convert beta-trigonometry to 
		coeff_correction = current.trc.return_coeff(angle_cur);
		//convert sz to numeric
		for (unsigned int i = 0; i < current.sz_power; i++)
			coeff_correction *= sz;
		//look for the same term
		it_map = (*it_vec).find(current.aop_c);
		if (it_map != it_vec->end()){
			it_map->second += coeff_correction*current.aop_c.coeff;
		}
		else {
			it_vec->insert({ current.aop_c, coeff_correction*current.aop_c.coeff });
		}

		it_vec++;//go to next point's storage
		angle_cur += angle_step;
	}
}

void AopXZRotateStorage::add_numerical_2_angles(AopXZRotate &current) {
	double angle_cur;
	double angle_cur_2 = angle_start_2;
	std::vector<std::vector<std::unordered_map<a_op_couple, double>>>::iterator it_vec_out = storage_numerical_2_angle.begin();
	std::vector<std::unordered_map<a_op_couple, double>>::iterator it_vec_in;
	std::unordered_map<a_op_couple, double>::iterator it_map;
	double coeff_correction;
	while (angle_cur_2 < angle_finish_2){
		it_vec_in = it_vec_out->begin();
		angle_cur = angle_start;
		while (angle_cur < angle_finish) {
			//convert beta-trigonometry to 
			coeff_correction = current.trc.return_coeff(angle_cur);
			coeff_correction *= current.trc2.return_coeff(angle_cur_2);
			//convert sz to numeric
			for (unsigned int i = 0; i < current.sz_power; i++)
				coeff_correction *= sz;
			//look for the same term
			it_map = (*it_vec_in).find(current.aop_c);
			if (it_map != it_vec_in->end()){
				it_map->second += coeff_correction*current.aop_c.coeff;
			}
			else {
				it_vec_in->insert({ current.aop_c, coeff_correction*current.aop_c.coeff });
			}
			it_vec_in++;//go to next point's storage
			angle_cur += angle_step;
		}
		it_vec_out++;
	}
}

void AopXZRotateStorage::add_operator(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos)
{
	if (operator_type == 0){ //select Sp'-term for  a-operator
		current.aop_c.names[new_operator_pos] = 'p';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1;	   //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Cos[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta/2]
			break;
		case 'p':
			current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		case 'm':
			current.aop_c.coeff *= -1;
			current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		}
	}
	else{ //select Sm'-term for a-operator
		current.aop_c.names[new_operator_pos] = 'm';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1;    //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
			break;
		case 'p':
			current.aop_c.coeff *= -1;
			current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		case 'm':
			current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		}
	}
}

void AopXZRotateStorage::add_operator_pi_zero(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos)
{
	if (operator_type == 0){ //select Sp'-term for  a-operator
		current.aop_c.names[new_operator_pos] = 'p';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1*Cos::getSignPiZero(t.nums[num]);	   //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Cos[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta/2]
			break;
		case 'p':
			current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		case 'm':
			current.aop_c.coeff *= -1;
			current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		}
	}
	else{ //select Sm'-term for a-operator
		current.aop_c.names[new_operator_pos] = 'm';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1*Cos::getSignPiZero(t.nums[num]);    //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
			break;
		case 'p':
			current.aop_c.coeff *= -1;
			current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		case 'm':
			current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		}
	}
}

void AopXZRotateStorage::add_operator_zero_pi_2_angles(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos)
{
	bool second_type = Cos::getSignZeroPi(t.nums[num]) == -1;
	if (operator_type == 0){ //select Sp'-term for  a-operator
		current.aop_c.names[new_operator_pos] = 'p';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1;	   //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.add_trig(2, 1, second_type); //Cos[beta/2]
			current.add_trig(3, 1, second_type); //Sin[beta/2]
			break;
		case 'p':
			current.add_trig(2, 2, second_type); //Cos^2[beta/2]
			break;
		case 'm':
			current.aop_c.coeff *= -1;
			current.add_trig(3, 2, second_type); //Sin^2[beta/2]
			break;
		}
	}
	else{ //select Sm'-term for a-operator
		current.aop_c.names[new_operator_pos] = 'm';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1;    //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.add_trig(2, 1, second_type); //Cos[beta/2]
			current.add_trig(3, 1, second_type); //Sin[beta/2]
			break;
		case 'p':
			current.aop_c.coeff *= -1;
			current.add_trig(3, 2, second_type);//Sin^2[beta/2]
			break;
		case 'm':
			current.add_trig(2, 2, second_type);//Cos^2[beta/2]
			break;
		}
	}
}

void AopXZRotateStorage::add_operator_antiferro(AopXZRotate &current, term &t,int num,int operator_type, int new_operator_pos)
{
	if (operator_type == 0){ //select Sp'-term for  a-operator
		current.aop_c.names[new_operator_pos] = 'p';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1;	   //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Cos[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta/2]
			break;
		case 'p':
			//current.aop.coeff *= 1;
			if (Cos::getSign(t.nums[num]) == 1)
				current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			else
				current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		case 'm':
			current.aop_c.coeff *= -1;
			if (Cos::getSign(t.nums[num]) == 1)
				current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			else
				current.trc.incCoeff(2, 2);//Cos^2[beta/2]				
			break;
		}
	}
	else{ //select Sm'-term for a-operator
		current.aop_c.names[new_operator_pos] = 'm';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -0.5;
			current.aop_c.coeff *= 2;    //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
			break;
		case 'p':
			current.aop_c.coeff *= -1;
			if (Cos::getSign(t.nums[num]) == 1)
				current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			else
				current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		case 'm':
			if (Cos::getSign(t.nums[num]) == 1)
				current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			else
				current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		}
	}
}

void AopXZRotateStorage::add_operator_2_sublattices(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos)
{
	if (operator_type == 0){ //select Sp'-term for  a-operator
		current.aop_c.names[new_operator_pos] = 'p';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1 * Cos::getSign(t.nums[num]);	   //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Cos[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta/2]
			break;
		case 'p':
			current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		case 'm':
			current.aop_c.coeff *= -1;
			current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		}
	}
	else{ //select Sm'-term for a-operator
		current.aop_c.names[new_operator_pos] = 'm';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1 * Cos::getSign(t.nums[num]);    //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Cos[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta/2]
			break;
		case 'p':
			current.aop_c.coeff *= -1;
			current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		case 'm':
			current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		}
	}
}

void AopXZRotateStorage::add_operator_ZY2_sublattices(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos)
{
	if (operator_type == 0){ //select Sp'-term for  a-operator
		current.aop_c.names[new_operator_pos] = 'p';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= -1 * Cos::getSign(t.nums[num]);	   //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Cos[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta/2]
			break;
		case 'p':
			current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		case 'm':
			current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		}
	}
	else{ //select Sm'-term for a-operator
		current.aop_c.names[new_operator_pos] = 'm';
		switch (t.ops[num])	{
		case 'z':
			current.aop_c.coeff *= Cos::getSign(t.nums[num]);    //-Sin[beta]/2=-1*Cos[beta/2]*Sin[beta/2]
			current.trc.incCoeff(2, 1);//Cos[beta/2]
			current.trc.incCoeff(3, 1);//Sin[beta/2]
			break;
		case 'p':
			current.trc.incCoeff(3, 2);//Sin^2[beta/2]
			break;
		case 'm':
			current.trc.incCoeff(2, 2);//Cos^2[beta/2]
			break;
		}
	}
}

void AopXZRotateStorage::ConvertToBilinear(term t) {
	if (t.len == -1) return; //no action in C-case
	AopXZRotate current;
	//select only z-terms
	current.aop_c.names[0] = 'p';
	current.aop_c.names[1] = 'm';
	current.aop_c.dx = 0;
	current.aop_c.dy = 0;
	current.aop_c.coeff = t.value;
	current.sz_power = t.len - 1; //amount of sz multipliers equal length-1

	//for correct sign applies multiplier "-2*sz". It's equal to sign of (a+)(a) in case 
	// sz=+-1/2;
	current.sz_power++;
	current.aop_c.coeff *= -2;
	//end sign 


	for (unsigned int i = 0; i < t.len; i++) {
		if (t.ops[i] != 'z')
		{
			//inc amount of Sin[beta]=2*Sin[beta/2]*Cos[beta/2] (p/m case)
			current.trc.incCoeff(2, 1);
			current.trc.incCoeff(3, 1);
			current.aop_c.coeff *= 2;
		}
		else
			current.trc.incCoeff(0, 1);//inc amount of Cos[beta] (z case)
	}

	if (analytical_mode == true) {
		std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
		if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
		else storage.insert({ current, current.aop_c.coeff }); //add new element;
	}
	else {
		add_numerical(current);
	}


	//check if enough operators fo pm case
	if (t.len < 2) return; //exit in case when not enough

	//select all combination of 2 operators as pm
	std::vector<std::pair<int, int>> pairs = generate_pairs(t.len);
	for (auto elem : pairs) {
		//Sz=Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
		//Sm=Sin[beta]*Sz'-Sin^2[beta/2]*Sp'+cos^2[beta/2]*Sm'

		//check all possible combinations 
		for (unsigned int i = 0; i < 4; i++) { //try all combinations 4: pp,pm,mp,mm
			current.clear();
			current.aop_c.coeff = t.value;

			//add shift between 2 nodes for K-type Cosine
			//output alwayas has positive dx
			Cos::findArbitraryCos(t.nums[elem.first], t.nums[elem.second], current.aop_c.dx, current.aop_c.dy);

			//first operator
			add_operator(current, t, elem.first, arr[i][0], 0);

			//second operator
			add_operator(current, t, elem.second, arr[i][1], 1);

			//add trig factors from "z-replaced" terms
			int skip = elem.first; //first element that not replaced to "z"
			for (unsigned int i = 0; i < t.len; i++) {
				if (i != skip) {
					switch (t.ops[i]) {
						//Sz=Cos[beta]*Sz'+...
						//Sp=Sin[beta]*Sz'+...
						//Sp=Sin[beta]*Sz'+...
					case 'z':
						current.trc.incCoeff(0, 1);
						break;
					case 'p':
						current.aop_c.coeff *= 2;	//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'm':
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					}
				}
				else
					skip = elem.second; //when find and skip first element, we change skip-variable to the number of second one
			}

			current.aop_c.coeff /= t.len; //according to duplication of translationally invariant routes
			current.sz_power = t.len - 2;

			current.aop_c.check();
			if (analytical_mode == true){
				std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
				if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
				else storage.insert({ current, current.aop_c.coeff }); //add new element;
			}
			else{
				add_numerical(current);
			}
		}
	}
}

void AopXZRotateStorage::ConvertToBilinearPiZero(term t) {
	if (t.len == -1) return; //no action in C-case
	AopXZRotate current;
	//select only z-terms
	current.aop_c.names[0] = 'p';
	current.aop_c.names[1] = 'm';
	current.aop_c.dx = 0;
	current.aop_c.dy = 0;
	current.aop_c.coeff = t.value;
	current.sz_power = t.len - 1; //amount of sz multipliers equal length-1

	//for correct sign applies multiplier "-2*sz". It's equal to sign of (a+)(a) in case 
	// sz=+-1/2;
	current.sz_power++;
	current.aop_c.coeff *= -2;
	//end sign 


	for (unsigned int i = 0; i < t.len; i++) {
		if (t.ops[i] != 'z')
		{
			//inc amount of Sin[beta]=2*Sin[beta/2]*Cos[beta/2] (p/m case)
			current.trc.incCoeff(2, 1);
			current.trc.incCoeff(3, 1);
			current.aop_c.coeff *= 2;
			//sign
			current.aop_c.coeff *= Cos::getSignPiZero(t.nums[i]);
			//end sign
		}
		else{
			current.trc.incCoeff(0, 1);//inc amount of Cos[beta] (z case)
		}
	}

	if (analytical_mode == true) {
		std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
		if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
		else storage.insert({ current, current.aop_c.coeff }); //add new element;
	}
	else {
		add_numerical(current);
	}


	//check if enough operators fo pm case
	if (t.len < 2) return; //exit in case when not enough

	//select all combination of 2 operators as pm
	std::vector<std::pair<int, int>> pairs = generate_pairs(t.len);
	for (auto elem : pairs) {
		//Exp(i*R*(0,pi))==1
		//Sz=Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
		//Sm=Sin[beta]*Sz'-Sin^2[beta/2]*Sp'+Cos^2[beta/2]*Sm'

		//Exp(i*R*(0,pi))==-1
		//Sz=Cos[beta]*Sz'+0.5*Sin[beta]*Sp'+0.5*Sin[beta]*Sm'
		//Sp=-Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
		//Sm=-Sin[beta]*Sz'-Sin^2[beta/2]*Sp'+Cos^2[beta/2]*Sm'

		//check all possible combinations 
		for (unsigned int i = 0; i < 4; i++) { //try all combinations 4: pp,pm,mp,mm
			current.clear();
			current.aop_c.coeff = t.value;

			//add shift between 2 nodes for K-type Cosine
			//output alwayas has positive dx
			Cos::findArbitraryCos(t.nums[elem.first], t.nums[elem.second], current.aop_c.dx, current.aop_c.dy);

			//first operator
			add_operator_pi_zero(current, t, elem.first, arr[i][0], 0);

			//second operator
			add_operator_pi_zero(current, t, elem.second, arr[i][1], 1);

			//add trig factors from "z-replaced" terms
			int skip = elem.first; //first element that not replaced to "z"
			for (unsigned int i = 0; i < t.len; i++) {
				if (i != skip) {
					switch (t.ops[i]) {
					case 'z':
						current.trc.incCoeff(0, 1);
						break;
					case 'p':
						current.aop_c.coeff *= Cos::getSignPiZero(t.nums[i]);
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'm':
						current.aop_c.coeff *= Cos::getSignPiZero(t.nums[i]);
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					}
				}
				else
					skip = elem.second; //when find and skip first element, we change skip-variable to the number of second one
			}

			current.aop_c.coeff /= t.len; //according to duplication of translationally invariant routes
			current.sz_power = t.len - 2;

			current.aop_c.check();
			if (analytical_mode == true){
				std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
				if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
				else storage.insert({ current, current.aop_c.coeff }); //add new element;
			}
			else{
				add_numerical(current);
			}
		}
	}
}

void AopXZRotateStorage::ConvertToBilinearZeroPi2Angles(term t) {
	if (t.len == -1) return; //no action in C-case
	

	AopXZRotate current;
	//select only z-terms
	current.aop_c.names[0] = 'p';
	current.aop_c.names[1] = 'm';
	current.aop_c.dx = 0;
	current.aop_c.dy = 0;
	current.aop_c.coeff = t.value;
	current.sz_power = t.len - 1; //amount of sz multipliers equal length-1

	//for correct sign applies multiplier "-2*sz". It's equal to sign of (a+)(a) in case 
	// sz=+-1/2;
	current.sz_power++;
	current.aop_c.coeff *= -2;
	//end sign 

	bool second_type;
	for (unsigned int i = 0; i < t.len; i++) {
		second_type = Cos::getSignZeroPi(t.nums[i]) == -1;
		if (t.ops[i] != 'z')
		{
			//inc amount of Sin[beta]=2*Sin[beta/2]*Cos[beta/2] (p/m case)
			current.add_trig(2, 1,second_type);
			current.add_trig(3, 1, second_type);
			current.aop_c.coeff *= 2;
		}
		else{
			current.add_trig(0, 1, second_type);//inc amount of Cos[beta] (z case)
		}
	}

	if (analytical_mode == true) {
		std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
		if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
		else storage.insert({ current, current.aop_c.coeff }); //add new element;
	}
	else {
		add_numerical(current);
	}


	//check if enough operators for pm case
	if (t.len < 2) return; //exit in case when not enough

	//select all combination of 2 operators as pm
	std::vector<std::pair<int, int>> pairs = generate_pairs(t.len);
	for (auto elem : pairs) {
		//Sz=Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
		//Sm=Sin[beta]*Sz'-Sin^2[beta/2]*Sp'+Cos^2[beta/2]*Sm'

		//check all possible combinations 
		for (unsigned int i = 0; i < 4; i++) { //try all combinations 4: pp,pm,mp,mm
			current.clear();
			current.aop_c.coeff = t.value;

			//add shift between 2 nodes for K-type Cosine
			//output alwayas has positive dx
			Cos::findArbitraryCos(t.nums[elem.first], t.nums[elem.second], current.aop_c.dx, current.aop_c.dy);

			//first operator
			add_operator_zero_pi_2_angles(current, t, elem.first, arr[i][0], 0);

			//second operator
			add_operator_zero_pi_2_angles(current, t, elem.second, arr[i][1], 1);

			//add trig factors from "z-replaced" terms
			int skip = elem.first; //first element that not replaced to "z"
			for (unsigned int i = 0; i < t.len; i++) {
				if (i != skip) {
					switch (t.ops[i]) {
					case 'z':
						current.add_trig(0, 1,second_type);
						break;
					case 'p':
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.add_trig(2, 1, second_type); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.add_trig(3, 1, second_type); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'm':
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.add_trig(2, 1, second_type); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.add_trig(3, 1, second_type); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					}
				}
				else
					skip = elem.second; //when find and skip first element, we change skip-variable to the number of second one
			}

			current.aop_c.coeff /= t.len; //according to duplication of translationally invariant routes
			current.sz_power = t.len - 2;

			current.aop_c.check();
			if (analytical_mode == true){
				std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
				if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
				else storage.insert({ current, current.aop_c.coeff }); //add new element;
			}
			else{
				add_numerical(current);
			}
		}
	}
}

void AopXZRotateStorage::ConvertToBilinearAntiFerro(term t) {
	if (t.len == -1) return; //no action in C-case
	AopXZRotate current;
	//select only z-terms
	current.aop_c.names[0] = 'p';
	current.aop_c.names[1] = 'm';
	current.aop_c.dx = 0;
	current.aop_c.dy = 0;
	current.aop_c.coeff = t.value;
	current.sz_power = t.len - 1; //amount of sz multipliers equal length-1

	//for correct sign applies multiplier "-2*sz". It's equal to sign of (a+)(a) in case 
	// sz=+-1/2;
	current.sz_power++;
	current.aop_c.coeff *= -2;
	//end sign 


	for (unsigned int i = 0; i < t.len; i++) {
		if (t.ops[i] != 'z')
		{
			//inc amount of Sin[beta]=2*Sin[beta/2]*Cos[beta/2] (p/m case)
			current.trc.incCoeff(2, 1);
			current.trc.incCoeff(3, 1);
			current.aop_c.coeff *= 2;
		}
		else{
			//sign
			current.aop_c.coeff *= Cos::getSign(t.nums[i]);
			//end sign
			current.trc.incCoeff(0, 1);//inc amount of Cos[beta] (z case)
		}			
	}

	if (analytical_mode == true) {
		std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
		if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
		else storage.insert({ current, current.aop_c.coeff }); //add new element;
	}
	else {
		add_numerical(current);
	}


	//check if enough operators fo pm case
	if (t.len < 2) return; //exit in case when not enough

	//select all combination of 2 operators as pm
	std::vector<std::pair<int, int>> pairs = generate_pairs(t.len);
	for (auto elem : pairs) {
		//Exp(i*R*k)==1
		//Sz=Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
		//Sp=Sin[beta]*Sz'-Sin^2[beta/2]*Sp'+Cos^2[beta/2]*Sm'
		
		//Exp(i*R*k)==-1
		//Sz=-Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Sin^2[beta/2]*Sp'+(-1)*Cos^2[beta/2]*Sm'
		//Sp=Sin[beta]*Sz'-Cos^2[beta/2]*Sp'+Sin^2[beta/2]*Sm'

		//check all possible combinations 
		for (unsigned int i = 0; i < 4; i++) { //try all combinations 4: pp,pm,mp,mm
			current.clear();
			current.aop_c.coeff = t.value;

			//add shift between 2 nodes for K-type Cosine
			//output alwayas has positive dx
			Cos::findArbitraryCos(t.nums[elem.first], t.nums[elem.second], current.aop_c.dx, current.aop_c.dy);

			//first operator
			add_operator_antiferro(current, t, elem.first, arr[i][0], 0);

			//second operator
			add_operator_antiferro(current, t, elem.second, arr[i][1], 1);
			
			//add trig factors from "z-replaced" terms
			int skip = elem.first; //first element that not replaced to "z"
			for (unsigned int i = 0; i < t.len; i++) {
				if (i != skip) {
					switch (t.ops[i]) {
					case 'z':
						current.trc.incCoeff(0, 1);
						current.aop_c.coeff *= Cos::getSign(t.nums[i]);
						break;
					case 'p':
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'm':
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					}
				}
				else
					skip = elem.second; //when find and skip first element, we change skip-variable to the number of second one
			}

			current.aop_c.coeff /= t.len; //according to duplication of translationally invariant routes
			current.sz_power = t.len - 2;

			current.aop_c.check();
			if (analytical_mode == true){
				std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
				if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
				else storage.insert({ current, current.aop_c.coeff }); //add new element;
			}
			else{
				add_numerical(current);
			}
		}
	}
}

void AopXZRotateStorage::ConvertToBilinear2sublattices(term t) {
	if (t.len == -1) return; //no action in C-case
	AopXZRotate current;
	//select only z-terms
	current.aop_c.names[0] = 'p';
	current.aop_c.names[1] = 'm';
	current.aop_c.dx = 0;
	current.aop_c.dy = 0;
	current.aop_c.coeff = t.value;
	current.sz_power = t.len - 1; //amount of sz multipliers equal length-1

	//for correct sign applies multiplier "-2*sz". It's equal to sign of (a+)(a) in case 
	// sz=+-1/2;
	current.sz_power++;
	current.aop_c.coeff *= -2;
	//end sign 


	for (unsigned int i = 0; i < t.len; i++) {
		if (t.ops[i] != 'z')
		{
			//inc amount of Sin[beta]=2*Sin[beta/2]*Cos[beta/2] (p/m case)
			current.trc.incCoeff(2, 1);
			current.trc.incCoeff(3, 1);
			current.aop_c.coeff *= 2;
			//sign
			current.aop_c.coeff *= Cos::getSign(t.nums[i]);
			//end sign
		}
		else{
			current.trc.incCoeff(0, 1);//inc amount of Cos[beta] (z case)
		}
	}

	if (analytical_mode == true) {
		std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
		if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
		else storage.insert({ current, current.aop_c.coeff }); //add new element;
	}
	else {
		add_numerical(current);
	}


	//check if enough operators fo pm case
	if (t.len < 2) return; //exit in case when not enough

	//select all combination of 2 operators as pm
	std::vector<std::pair<int, int>> pairs = generate_pairs(t.len);
	for (auto elem : pairs) {
		//Exp(i*R*k)==1
		//Sz=Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
		//Sp=Sin[beta]*Sz'-Sin^2[beta/2]*Sp'+Cos^2[beta/2]*Sm'

		//Exp(i*R*k)==-1
		//Sz=-Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Sin^2[beta/2]*Sp'+(-1)*Cos^2[beta/2]*Sm'
		//Sp=Sin[beta]*Sz'-Cos^2[beta/2]*Sp'+Sin^2[beta/2]*Sm'

		//check all possible combinations 
		for (unsigned int i = 0; i < 4; i++) { //try all combinations 4: pp,pm,mp,mm
			current.clear();
			current.aop_c.coeff = t.value;

			//add shift between 2 nodes for K-type Cosine
			//output alwayas has positive dx
			Cos::findArbitraryCos(t.nums[elem.first], t.nums[elem.second], current.aop_c.dx, current.aop_c.dy);

			//first operator
			add_operator_2_sublattices(current, t, elem.first, arr[i][0], 0);

			//second operator
			add_operator_2_sublattices(current, t, elem.second, arr[i][1], 1);

			//add trig factors from "z-replaced" terms
			int skip = elem.first; //first element that not replaced to "z"
			for (unsigned int i = 0; i < t.len; i++) {
				if (i != skip) {
					switch (t.ops[i]) {
					case 'z':
						current.trc.incCoeff(0, 1);
						break;
					case 'p':
						current.aop_c.coeff *= Cos::getSign(t.nums[i]);
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'm':
						current.aop_c.coeff *= Cos::getSign(t.nums[i]);
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					}
				}
				else
					skip = elem.second; //when find and skip first element, we change skip-variable to the number of second one
			}

			current.aop_c.coeff /= t.len; //according to duplication of translationally invariant routes
			current.sz_power = t.len - 2;

			current.aop_c.check();
			if (analytical_mode == true){
				std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
				if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
				else storage.insert({ current, current.aop_c.coeff }); //add new element;
			}
			else{
				add_numerical(current);
			}
		}
	}
}

void AopXZRotateStorage::ConvertToBilinearZY2sublattices(term t) {
	if (t.len == -1) return; //no action in C-case
	AopXZRotate current;
	//select only z-terms
	current.aop_c.names[0] = 'p';
	current.aop_c.names[1] = 'm';
	current.aop_c.dx = 0;
	current.aop_c.dy = 0;
	current.aop_c.coeff = t.value;
	current.sz_power = t.len - 1; //amount of sz multipliers equal length-1

	//for correct sign applies multiplier "-2*sz". It's equal to sign of (a+)(a) in case 
	// sz=+-1/2;
	current.sz_power++;
	current.aop_c.coeff *= -2;
	//end sign 


	for (unsigned int i = 0; i < t.len; i++) {
		if (t.ops[i] == 'p')
		{
			//inc amount of Sin[beta]=2*Sin[beta/2]*Cos[beta/2] (p/m case)
			current.trc.incCoeff(2, 1);
			current.trc.incCoeff(3, 1);
			current.aop_c.coeff *= 2;
			//sign
			current.aop_c.coeff *= Cos::getSign(t.nums[i]);
			//end sign
			current.aop_c.i_power++;
			if (current.aop_c.i_power == 2){
				current.aop_c.i_power = 0;
				current.aop_c.coeff *= -1;
			}

		}
		else if (t.ops[i] == 'm')
		{
			//inc amount of Sin[beta]=2*Sin[beta/2]*Cos[beta/2] (p/m case)
			current.trc.incCoeff(2, 1);
			current.trc.incCoeff(3, 1);
			current.aop_c.coeff *= 2;
			//sign
			current.aop_c.coeff *= -1 * Cos::getSign(t.nums[i]);
			//end sign
			current.aop_c.i_power++;
			if (current.aop_c.i_power == 2){
				current.aop_c.i_power = 0;
				current.aop_c.coeff *= -1;
			}
		}
		else{
			current.trc.incCoeff(0, 1);//inc amount of Cos[beta] (z case)
		}
	}

	if (analytical_mode == true) {
		std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
		if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
		else storage.insert({ current, current.aop_c.coeff }); //add new element;
	}
	else {
		add_numerical(current);
	}


	//check if enough operators fo pm case
	if (t.len < 2) return; //exit in case when not enough

	//select all combination of 2 operators as pm
	std::vector<std::pair<int, int>> pairs = generate_pairs(t.len);
	for (auto elem : pairs) {
		//Exp(i*R*k)==1
		//Sz=Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
		//Sp=Sin[beta]*Sz'-Sin^2[beta/2]*Sp'+Cos^2[beta/2]*Sm'

		//Exp(i*R*k)==-1
		//Sz=-Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
		//Sp=Sin[beta]*Sz'+Sin^2[beta/2]*Sp'+(-1)*Cos^2[beta/2]*Sm'
		//Sp=Sin[beta]*Sz'-Cos^2[beta/2]*Sp'+Sin^2[beta/2]*Sm'

		//check all possible combinations 
		for (unsigned int i = 0; i < 4; i++) { //try all combinations 4: pp,pm,mp,mm
			current.clear();
			current.aop_c.i_power = 0;
			
			current.aop_c.coeff = t.value;

			//add shift between 2 nodes for K-type Cosine
			//output alwayas has positive dx
			Cos::findArbitraryCos(t.nums[elem.first], t.nums[elem.second], current.aop_c.dx, current.aop_c.dy);

			//first operator
			add_operator_2_sublattices(current, t, elem.first, arr[i][0], 0);

			//second operator
			add_operator_2_sublattices(current, t, elem.second, arr[i][1], 1);

			//add trig factors from "z-replaced" terms
			int skip = elem.first; //first element that not replaced to "z"
			for (unsigned int i = 0; i < t.len; i++) {
				if (i != skip) {
					switch (t.ops[i]) {
					case 'z':
						current.trc.incCoeff(0, 1);
						break;
					case 'p':
						current.aop_c.coeff *= Cos::getSign(t.nums[i]);
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.aop_c.i_power++;
						if (current.aop_c.i_power == 2){
							current.aop_c.i_power = 0;
							current.aop_c.coeff *= -1;
						}
						break;
					case 'm':
						current.aop_c.coeff *= -1 * Cos::getSign(t.nums[i]);
						current.aop_c.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.aop_c.i_power++;
						if (current.aop_c.i_power == 2){
							current.aop_c.i_power = 0;
							current.aop_c.coeff *= -1;
						}
						break;
					}
				}
				else
					skip = elem.second; //when find and skip first element, we change skip-variable to the number of second one
			}

			current.aop_c.coeff /= t.len; //according to duplication of translationally invariant routes
			current.sz_power = t.len - 2;

			current.aop_c.check();
			if (analytical_mode == true){
				std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
				if (it != storage.end()) it->second += current.aop_c.coeff; //add value to existed element
				else storage.insert({ current, current.aop_c.coeff }); //add new element;
			}
			else{
				add_numerical(current);
			}
		}
	}
}

void AopXZRotateStorage::print(std::ofstream &F) {
	F.precision(8);
	F << std::fixed;
	for (auto elem = storage.begin(); elem != storage.end(); ++elem)
	{
		if (abs(elem->second) < 0.00000001) return;
		F << "+" << elem->second << "*";
		elem->first.aop_c.printAterm(F, matrix, matrix_size, false);
		elem->first.trc.printTrigCoeff(F);
		switch (elem->first.sz_power){
		case 0: break;
		case 1: F << "*sz"; break;
		default: F << "*(sz^" << elem->first.sz_power << ")";
		}

	}
}

void AopXZRotateStorage::print_numerical(std::vector<std::ofstream*> &Fs) {
	std::vector<std::ofstream*>::iterator it_of = Fs.begin();
	std::vector<std::unordered_map<a_op_couple, double>>::iterator it_map = storage_numerical.begin();
	//test declaration
	int streamnum = -1,elnum=-1;
	//end test declaration
	while (it_of < Fs.end()){
		streamnum++; //test
		elnum = -1;
		for (auto &elem : (*it_map)){
			elnum++; //test

			//if (elem.first.i_power==0)	
			//	(*(*it_of)) << "+" << elem.second << "*";
			//else if (abs(elem.second) > 0.00000001)
			//	std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAA \n AAAAAAAAAAAAAAAAAAAAAAAAA\n AAAAAAAAAAAAAAAAAAAA\n!!!!!!!!!!!!!!!!";
			if (abs(elem.second) > 0.000000001){
				(*(*it_of)) << "+" << elem.second << "*";
				elem.first.printAterm((*(*it_of)), matrix, matrix_size, false);
			}
		}
		it_of++;
		it_map++;
	}
}

void AopXZRotateStorage::print_numerical_2_angles(std::vector<std::vector<std::ofstream*>> &Fs) {
	std::vector<std::vector<std::ofstream*>>::iterator it_ofstr_out = Fs.begin();
	std::vector<std::ofstream*>::iterator it_ofstr_in;

	std::vector<std::vector<std::unordered_map<a_op_couple, double>>>::iterator it_vec_out = storage_numerical_2_angle.begin();
	std::vector<std::unordered_map<a_op_couple, double>>::iterator it_vec_in;
	//test declaration
	int streamnum = -1, elnum = -1;
	//end test declaration
	while (it_ofstr_out < Fs.end()){
		it_vec_in = it_vec_out->begin();
		it_ofstr_in = it_ofstr_out->begin();
		while (it_ofstr_in < it_ofstr_out->end()){
			//streamnum++; //test
			//elnum = -1;
			for (auto &elem : (*it_vec_in)){
				//elnum++; //test
				//if (elem.first.i_power==0)	
				//	(*(*it_of)) << "+" << elem.second << "*";
				//else if (abs(elem.second) > 0.00000001)
				//	std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAA \n AAAAAAAAAAAAAAAAAAAAAAAAA\n AAAAAAAAAAAAAAAAAAAA\n!!!!!!!!!!!!!!!!";
				if (abs(elem.second) > 0.000000001){
					(**it_ofstr_in) << "+" << elem.second << "*";
					elem.first.printAterm(**it_ofstr_in, matrix, matrix_size, false);
				}
			}
			it_ofstr_in++;
			it_vec_in++;
		}
		it_vec_out++;
		it_ofstr_out++;
	}
}