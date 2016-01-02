// 9 Collect Third Terms.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


#include "converter.h"
#include "term.h"
#include "a_op.h"
#include "ground_energy.h"
#include "trig_coefficients.h"
#include "a_op_xz_rotate.h"

using namespace std;


const int min_op_amount = 4;
const string delim = "\\";
const string config_dir = "config";
const string inp_res = "input_rot";
const string out_res = "Results_rot";
const string out_file_end = "_rot.txt";

//Trigonometry coefficients for XZ spin-rotation



int arr[4][2] = { { 0, 0 }, { 0, 1 }, { 1, 0 }, { 1, 1 } };



//Hash functions for unordered maps
namespace std{
	template <>
	struct hash<a_op>
	{
		std::size_t operator()(const a_op& k) const
		{
			std::ostringstream out;
			out << " " << k.coeff << " " << k.n << " " << k.coeff << " ";
			for (auto elem : k.names)
				out << elem << " ";
			for (auto elem : k.node)
				out << elem << " ";
			std::string s = out.str();
			return (std::hash<string>()(s));
		}
	};
}
namespace std{
	template <>
	struct hash<AopXZRotate>
	{
		std::size_t operator()(const AopXZRotate& k) const
		{
			std::ostringstream out;
			out << " " << k.aop.coeff << " " << k.aop.n << " " << k.aop.coeff << " ";
			for (auto elem : k.aop.names)
				out << elem << " ";
			for (auto elem : k.aop.node)
				out << elem << " ";
			out << k.sz_power << " ";

			out << k.trc.getPowers();
			std::string s = out.str();

			return (std::hash<string>()(s));
		}
	};
}

//storage class that contain all current elements

class AopXZRotateStorage{
public: 
	//members
	std::unordered_map<AopXZRotate, double> storage; //for terms with trigonometry coeffs
	std::vector<std::unordered_map<a_op, double>> storage_numerical;//for terms with numeric values only
	int **matrix; //matrix vith node numbers
	int matrix_size;
	bool analytical_mode; //true- analytical, false -numerical;
	double angle_start;
	double angle_finish;
	double angle_step;
	double sz;

	//methods
	void set(int **Matrix, int Matrix_size, bool Mode, double Angle_start, double Angle_finish, double Angle_step,double Sz=-0.5)
	{
		matrix = Matrix;
		matrix_size = Matrix_size;
		analytical_mode = Mode;

		if (analytical_mode == false) {
			
			angle_start = Angle_start;
			angle_step = Angle_step;
			angle_finish = Angle_finish;
			double angle_cur = angle_start;
			while (angle_cur <= angle_finish) {
				storage_numerical.push_back(std::unordered_map<a_op, double>());
					angle_cur += Angle_step;
			}
		}
		sz = Sz;
	}
	
	void clearTerms()
	{
		storage.clear();
		for (unsigned int i = 0; i < storage_numerical.size();i++) {
			storage_numerical[i].clear();
		}
		storage_numerical.clear();
	}

	std::vector<std::pair<int, int>> generate_pairs(int n) {
		//Генерирует все возможные пары для данного терма
		std::vector<std::pair<int, int>> pairs;
		for (int i = 0; i < n; i++)	{
			for (int j = i + 1; j < n; j++) {
				pairs.push_back(std::pair<int, int>::pair(i, j));
			}
		}
		return pairs;
	}

	void add_numerical(AopXZRotate &current) {
		double angle_cur = angle_start;
		std::vector<std::unordered_map<a_op, double>>::iterator it_vec = storage_numerical.begin();
		std::unordered_map<a_op, double>::iterator it_map;
		double coeff_correction;
		while (angle_cur < angle_finish) {
			//convert beta-trigonometry to 
			coeff_correction = current.trc.return_coeff(angle_cur);
			//convert sz to numeric
			for (unsigned int i = 0; i < current.sz_power; i++)
				coeff_correction *= sz;
			//look for the same term
			(*it_vec).find(current.aop);
			it_map = (*it_vec).find(current.aop);
			if (it_map != it_vec->end()){
				it_map->second += coeff_correction*current.aop.coeff;
			}
			else {
				it_vec->insert({ current.aop, coeff_correction*current.aop.coeff });
			}

			it_vec++;//go to next point's storage
			angle_cur += angle_step;
		}
	}

	void ConvertToBilinear(term t) {
		if (t.len == -1) return; //no action in C-case
		AopXZRotate current;
		//select only z-terms
		current.aop.n = 2;
		current.aop.names.push_back('p');
		current.aop.names.push_back('m');
		current.aop.node.push_back(0);
		current.aop.node.push_back(0);
		current.aop.coeff = t.value;
		current.sz_power = t.len - 1; //amount of sz multipliers equal length-1

		//for correct sign applies multiplier "-2*sz". It's equal to sign of (a+)(a) in case 
		// sz=+-1/2;
		current.sz_power++;
		current.aop.coeff *= -2;
		//end sign 

		
		for (unsigned int i = 0; i < t.len; i++) {
			if (t.ops[i] != 'z')
			{
				//inc amount of Sin[beta]=2*Sin[beta/2]*Cos[beta/2] (p/m case)
				current.trc.incCoeff(2, 1);
				current.trc.incCoeff(3, 1); 
				current.aop.coeff *= 2;
			}
			else
				current.trc.incCoeff(0,1);//inc amount of Cos[beta] (z case)
		}
		if (analytical_mode == true) {
			std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
			if (it != storage.end()) it->second += current.aop.coeff; //add value to existed element
			else storage.insert({ current, current.aop.coeff }); //add new element;
		}
		else {
			add_numerical(current);
		}
		
		
		//check if enough operators fo pm case
		if (t.len < 2) return; //exit in case when not enough

		//select all combination of 2 operators as pm
		std::vector<std::pair<int, int>> pairs = generate_pairs(t.len);
		for (auto elem:pairs) {
			//Sz=Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
			//Sp=Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
			//Sp=Sin[beta]*Sz'-Cos^2[beta/2]*Sp'-Sin^2[beta/2]*Sm'
			
			//check all possible combinations 
			for (unsigned int i = 0; i < 4; i++) { //try all combinations 4: pp,pm,mp,mm
				current.clear();
				current.aop.n = 2;
				current.aop.coeff = t.value;
				//first operator
				current.aop.node.push_back(t.nums[elem.first]); //set first node-number
				if (arr[i][0] == 0){ //select Sp' -term for first a-operator
					current.aop.names.push_back('p');
					switch (t.ops[elem.first])	{
					case 'z':
						current.aop.coeff *= -0.5;
						current.aop.coeff *= 2;	   //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'p':
						//current.aop.coeff *= 1;
						current.trc.incCoeff(2,2);//Cos^2[beta/2]
						break;
					case 'm':
						current.aop.coeff *= -1;
						current.trc.incCoeff(3,2);//Sin^2[beta/2]
						break;
					}
				}
				else{ //select Sm-term for first a-operator
					current.aop.names.push_back('m');
					switch (t.ops[elem.first])	{
					case 'z':
						current.aop.coeff *= -0.5;
						current.aop.coeff *= 2;    //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'p':
						current.aop.coeff *= -1;
						current.trc.incCoeff(3,2);//Sin^2[beta/2]
						break;
					case 'm':
						//current.aop.coeff *= 1;
						current.trc.incCoeff(2,2);//Cos^2[beta/2]
						break;
					}
				}

				//second operator
				current.aop.node.push_back(t.nums[elem.second]); //set first node-number

				if (arr[i][1] == 0){ //select Sp -term for second a-operator
					current.aop.names.push_back('p');
					switch (t.ops[elem.second])	{
					case 'z':
						current.aop.coeff *= -0.5;
						current.aop.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'p':
						//current.aop.coeff *= 1;
						current.trc.incCoeff(2,2);//Cos^2[beta/2]
						break;
					case 'm':
						current.aop.coeff *= -1;
						current.trc.incCoeff(3,2);//Sin^2[beta/2]
						break;
					}
				}
				else{ //select Sm-term for second a-operator
					current.aop.names.push_back('m');
					switch (t.ops[elem.second])	{
					case 'z':
						current.aop.coeff *= -0.5;
						current.aop.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
						break;
					case 'p':
						current.aop.coeff *= -1;
						current.trc.incCoeff(3,2);//Sin^2[beta/2]
						break;
					case 'm':
						//current.aop.coeff *= 1;
						current.trc.incCoeff(2,2);//Cos^2[beta/2]
						break;
					}
				}
				//TODO add trig factors from "z-replaced" terms
				int skip = elem.first;
				for (unsigned int i = 0; i < t.len; i++) {
					if (i != elem.first) {
						switch (t.ops[i]) {
							//Sz=Cos[beta]*Sz'+(-0.5)*Sin[beta]*Sp'+(-0.5)*Sin[beta]*Sm'
							//Sp=Sin[beta]*Sz'+Cos^2[beta/2]*Sp'+(-1)*Sin^2[beta/2]*Sm'
							//Sp=Sin[beta]*Sz'-Cos^2[beta/2]*Sp'-Sin^2[beta/2]*Sm'
						case 'z': 
							current.trc.incCoeff(0, 1); 
							break;
						case 'p': 
							current.aop.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
							current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
							current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
							break;
						case 'm': 
							current.aop.coeff *= 2;		//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
							current.trc.incCoeff(2, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
							current.trc.incCoeff(3, 1); //Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
							break;
						}
					}
					else
						skip = elem.second;
				}

				current.aop.coeff /= t.len; //according to duplication of translationally invariant routes
				current.sz_power = t.len - 2;
				
				if (analytical_mode == true){
					std::unordered_map<AopXZRotate, double>::iterator it = storage.find(current);
					if (it != storage.end()) it->second += current.aop.coeff; //add value to existed element
					else storage.insert({ current, current.aop.coeff }); //add new element;
				}
				else{
					add_numerical(current);
				}
			}
		}
	}

	void print(ofstream &F)
	{
		F.precision(8);
		F << std::fixed;
		for (auto elem = storage.begin(); elem != storage.end();++elem)
		{
			if (abs(elem->second) < 0.00000001) return;
			F <<"+"<< elem->second << "*";
			elem->first.aop.printAterm(F, matrix, matrix_size,false);
			elem->first.trc.printTrigCoeff(F);
			switch (elem->first.sz_power){
			case 0: break;
			case 1: F << "*sz"; break;
			default: F << "*(sz^" << elem->first.sz_power << ")";
			}
			
		}
	}

	void print_numerical(std::vector<std::ofstream*> &Fs)
	{
		std::vector<std::ofstream*>::iterator it_of = Fs.begin();
		std::vector<std::unordered_map<a_op, double>>::iterator it_map = storage_numerical.begin();
		while (it_of < Fs.end()){
			for (auto &elem:(*it_map)){ 
				(*(*it_of)) << "+" << elem.second << "*";
				elem.first.printAterm((*(*it_of)), matrix, matrix_size, false);
			}
			it_of++;
			it_map++;
		}
	}
};

int** matrix;

void fillMatrix(int **matrix, int NNN)
{
	int dir = 0;
	int d = max(1 + ((NNN) / 2) * 2, 1 + (NNN - 2) * 2);
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
	//!------------Matrix output--------------------
	for (int i = 0; i<d; i++)
	{
		for (int j = 0; j<d; j++)
		{
			cout << setw(3) << matrix[i][j] << " ";
		}
		cout << "\n";
	}
	//!-----------Output end
	ofstream nodenums("nodenums.txt", ios::out);

	nodenums << "{";
	for (int i = 0; i<d; i++)
	{
		nodenums << "{";
		for (int j = 0; j<d; j++)
		{
			if (j == d - 1)
				nodenums << matrix[i][j] << "}";
			else
				nodenums << matrix[i][j] << ",";
		}
		if (i<d - 1)
			nodenums << ",";
	}
	nodenums << "}";
	nodenums.close();


}

int _tmain(int argc, _TCHAR* argv[])
{

	cout << min_op_amount << "\n";
	vector<string> points;
	std::ostringstream fname;
	fname << config_dir << delim << "config.txt";
	ifstream config(fname.str(), ios::in);
	int num_points, min_order, max_order;
	string tmp;
	int a_amount,mode,analytical_mode;
	double angle_start, angle_step, angle_finish;
	getline(config, tmp);
	config >> num_points;
	getline(config, tmp);
	getline(config, tmp);
	config >> min_order >> max_order;
	getline(config, tmp);
	getline(config, tmp);
	config >> a_amount;
	getline(config, tmp);
	getline(config, tmp);
	config >> mode;//
	getline(config, tmp);
	getline(config, tmp);
	config >> analytical_mode;//
	getline(config, tmp);
	getline(config, tmp);
	config >> angle_start>> angle_step>> angle_finish;
	config.close();
	///init matrix
	int size = max(1 + ((max_order) / 2) * 2, 1 + (max_order - 2) * 2);
	int** m;
	m = new int*[size];
	for (int i = 0; i < size; i++)
		m[i] = new int[size];
	fillMatrix(m, max_order);
	//end init

	fname.str("");
	fname << config_dir << delim << "points.txt";
	ifstream fpoints(fname.str(), ios::in);
	string s, tmp_s;
	for (int i = 0; i < num_points; i++)
	{
		getline(fpoints, s);
		points.push_back(s);
	}
	fpoints.close();



	double temp_val;

	for (int i = 0; i < num_points; i++)
	{
		cout << i << "\n";
		double factor;
		///test
		double factor_edge = 0.56;
		cout << "TEST FACTOR-EDGE";
		//end test
		if (stof(points[i]) - factor_edge > 0) ///
			factor = 0.5;
		else
			factor = -0.5;
		converter conv1;
		groundEnergy ge;
		AopXZRotateStorage rotate_excitation_storage;
		term t1;

		ofstream out, out_energy, out_rot, out_rot_excitation;
		ofstream *out_rot_cur;

		fname.str("");
		fname << out_res << delim << "res_" << points[i] << "_" << max_order << "_" << a_amount << ".txt";
		out.open(fname.str(), ios::out);
		out << "{";
		fname.str("");
		fname << out_res << delim << "energy_res_" << points[i] << "_" << max_order << "_" << a_amount << ".txt";
		out_energy.open(fname.str(), ios::out);
		out_energy << "{";

		//init of rotate energy file
		fname.str("");
		fname << out_res << delim << "energy_rot_" << points[i] << "_" << max_order << "_" << a_amount << ".txt";
		out_rot.open(fname.str(), ios::out);
		out_rot << "{";
		//end init

		//init of rotate excitations file
		fname.str("");
		fname << out_res << delim << "excitations_rot_" << points[i] << "_" << max_order << "_" << a_amount << ".txt";
		out_rot_excitation.open(fname.str(), ios::out);
		out_rot_excitation << "{";
		//end init

		//init of array for numerical excitations files
		std::vector<std::ofstream*> out_rot_excitation_numericals;
		double angle_cur = angle_start;
		while (angle_cur<angle_finish) {
			fname.str("");
			fname << out_res << delim << points[i] << delim << "excitations_rot_" << points[i] << "_" << max_order << "_" << a_amount << "_" << angle_cur << ".txt";
			out_rot_cur = new ofstream();
			(*out_rot_cur).open(fname.str(), ios::out);
			(*out_rot_cur) << "{";
			(*out_rot_cur) << std::fixed;
			(*out_rot_cur).precision(8);
			out_rot_excitation_numericals.push_back(out_rot_cur);
			angle_cur += angle_step;
		}
		//end init

		conv1.set(factor, a_amount, m, size);
		
		ge.set(factor);
		for (int j = min_order; j <= max_order; j++)
		{
			cout << "Order " << j << "\n";
			conv1.clearTerms();
			ge.clearTerms();
			rotate_excitation_storage.clearTerms();
			//Test  
			cout << "\n\n CHANGE MIN OP TO 1!!!!!!!!!!!!!!\n\n";
			//end Test block
			for (int k = min_op_amount; k <= j; k++)
			{
				cout << "SubOrder: " << k << "\n";
				fname.str("");
				fname << "d:\\Andrew\\Practice\\!!!_Last Set\\6.2 Collect Terms\\6.2 Collect Terms" << delim << inp_res << delim << points[i] << "\\" << j << "_results_" << points[i] << "_" << k << out_file_end;

				ifstream cur(fname.str(), ios::in);

				t1.setOrder(j);
				rotate_excitation_storage.set(m, size, analytical_mode != 0, angle_start, angle_finish, angle_step);
				int str_amount = 0;
				while (!cur.eof())
				{
					str_amount++;
					if (str_amount % 100 == 0)
						cout << str_amount << " ";
					getline(cur, s);
					if (s.length() > 0)
					{
						istringstream iss;
						iss.str(s);
						iss >> tmp_s >> temp_val;


						t1.decompose(tmp_s, temp_val);
						switch (mode){
						case 0:
							if (t1.type)//for energy
								ge.addTerm(j, t1);
							break;
						case 1:
							ge.addTermRotation(j, t1);
							break;
						case 2:
							if (t1.len > 0)//for excitations
								conv1.decomposeTerm(t1);
							break;
						case 3:
							if (t1.len > 0)//for rotate excitations
								rotate_excitation_storage.ConvertToBilinear(t1);
							break;
						}
					}
				}
				cur.close();
			}


			std::vector<std::ofstream*>::iterator it = out_rot_excitation_numericals.begin();
			switch (mode){

			case 0:
				if (a_amount = 2)
				{
					out_energy << ge.returnE(j);
					if (j != max_order)
					{
						out << ",";
						out_energy << ",";
					}
				}
				break;
			case 1:
				ge.printTermRotation(out_rot, j);
				if (j != max_order)
					out_rot << ",";
				break;
			case 2:
				if (a_amount = 2)
				{
					conv1.convertToAop();
					if (!conv1.PrintAllAop(out))
						out << "0";
				}
				break;
			case 3:
				switch (analytical_mode){
				case 0:
					rotate_excitation_storage.print_numerical(out_rot_excitation_numericals);
					
					if (j != max_order){
						while (it < out_rot_excitation_numericals.end()) {
							(**it)<<",";
							it++;
						}
					}
						
					break;
				case 1:
					rotate_excitation_storage.print(out_rot_excitation);
					if (j != max_order)
						out_rot_excitation << ",";
					break;
				}
				break;
			}
			std::cout << "Done!\n";
		}
		out << "}";
		out_energy << "}";
		out_rot << "}";
		out_rot_excitation << "}";

		out.close();
		out_energy.close();
		out_rot.close();
		out_rot_excitation.close();
		std::vector<std::ofstream*>::iterator it = out_rot_excitation_numericals.begin();
		while (it < out_rot_excitation_numericals.end()) {
			(**it).close();
			delete (*it);
			it++;
		}
	}

	return 0;
}

