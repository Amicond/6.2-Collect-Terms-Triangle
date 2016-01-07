// 9 Collect Third Terms.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


#include "converter.h"
#include "term.h"
#include "a_op.h"
#include "a_op_couple.h"
#include "ground_energy.h"
#include "trig_coefficients.h"
#include "a_op_xz_rotate.h"
#include "a_op_xz_rotate_storage.h"

using namespace std;


const int min_op_amount = 1;
const string delim = "\\";
const string config_dir = "config";
const string inp_res = "input_rot";
const string out_res = "Results_rot";
const string out_file_end = "_rot.txt";


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
			return (std::hash<string>()(out.str()));
		}
	};
}


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

	/*
	int iii;//cosine test
	std::cout << cos(3.1415 / 3) << std::endl;
	std::cin >> iii;
	return 0;*/

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
	Cos::set(m, size);
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

		ofstream out, out_energy, out_rot, out_energy_rot_antiferro,out_rot_excitation;
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

		//init of rotate energy file antiferromagnet case
		fname.str("");
		fname << out_res << delim << "energy_rot_antiferro_" << points[i] << "_" << max_order << "_" << a_amount << ".txt";
		out_energy_rot_antiferro.open(fname.str(), ios::out);
		out_energy_rot_antiferro.precision(10);
		out_energy_rot_antiferro << fixed;
		out_energy_rot_antiferro << "{";
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
		
		if (mode==4) //antiferromagnet case
			ge.set(-0.5,true);
		else
			ge.set(factor, true); //all other cases
		for (int j = min_order; j <= max_order; j++)
		{
			cout << "Order " << j << "\n";
			conv1.clearTerms();
			ge.clearTerms();
			rotate_excitation_storage.clearTerms();
			//Test  
			//cout << "\n\n CHANGE MIN OP TO 1!!!!!!!!!!!!!!\n\n";
			//end Test block
			rotate_excitation_storage.set(m, size, analytical_mode != 0, angle_start, angle_finish, angle_step);
			for (int k = min_op_amount; k <= j; k++)
			{
				cout << "SubOrder: " << k << "\n";
				fname.str("");
				fname << "d:\\Andrew\\Practice\\!!!_Last Set\\6.2 Collect Terms\\6.2 Collect Terms" << delim << inp_res << delim << points[i] << "\\" << j << "_results_" << points[i] << "_" << k << out_file_end;

				ifstream cur(fname.str(), ios::in);

				t1.setOrder(j);
				
				int str_amount = 0;
				while (!cur.eof())
				{
					str_amount++;
					if (str_amount % 200000 == 0)
						cout << "Or:" << j << " Sub:" << k << "\n";
					if (str_amount % 5000 == 0)
						cout<<str_amount <<" ";
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
						case 4:
								ge.addTermRotationAntiferromagnet(j, t1);
							break;
						}
					}
				}
				cur.close();
			}

			//defenition is here, because it's impossible to determine iterator in the switch
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
			case 4:
				ge.printTermRotation(out_energy_rot_antiferro, j);
				if (j != max_order)
					out_energy_rot_antiferro << ",";
				break;
			}
			std::cout << "Done!\n";
		}
		out << "}";
		out.close();

		out_energy << "}";
		out_energy.close();

		out_rot << "}";
		out_rot.close();

		out_energy_rot_antiferro << "}";
		out_energy_rot_antiferro.close();

		out_rot_excitation << "}";
		out_rot_excitation.close();
		
		std::vector<std::ofstream*>::iterator it = out_rot_excitation_numericals.begin();
		while (it < out_rot_excitation_numericals.end()) {
			(**it) << "}";
			(**it).close();
			delete (*it);
			it++;
		}
	}

	return 0;
}

