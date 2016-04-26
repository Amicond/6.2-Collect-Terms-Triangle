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
#include "term_storage.h"

using namespace std;


const int min_op_amount = 1;
const string delim = "\\";
const string config_dir = "config";
const string inp_res = "input_rot";
const string out_res = "Results_rot";
const string out_file_end = "_rot.txt";
const string add_type = "zero_pi_2_types"; //"antiferro";


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


//mode=zero
//energy
class eval_results
{

public:
	std::ostringstream fname;
	std::vector<string> *points;

	
	int min_order, max_order;
	int** m;
	int size;

	double angle_start, angle_step, angle_finish;
	double angle_start2, angle_step2, angle_finish2;
	double factor;


	eval_results(vector<string> &Points, int** M,int Size)
	{
		points = &Points;
		max_order = 1;
		m = M;
		size = Size;
	}

	

	//mode==1 ?
	void eval_energy_rotate(int i)
	{
		groundEnergy ge;
		term t1;
		string s, tmp_s;
		double temp_val;
		ge.set(factor, true); //all other cases
		std::ofstream out_rot;
		fname.str("");
		fname << out_res << delim << "energy_rot_" << (*points)[i] << "_" << max_order << ".txt";
		out_rot.open(fname.str(), ios::out);
		out_rot.precision(10);
		out_rot << fixed;
		out_rot << "{";
		for (int j = min_order; j <= max_order; j++)
		{
			cout << "Order " << j << "\n";
			ge.clearTerms();
			for (int k = min_op_amount; k <= j; k++)
			{
				cout << "SubOrder: " << k << "\n";
				fname.str("");
				fname << "d:\\Andrew\\Practice\\!!!_Last Set\\6.2 Collect Terms\\6.2 Collect Terms" << delim << inp_res << delim << (*points)[i] << "\\" << j << "_results_" << (*points)[i] << "_" << k << out_file_end;
				ifstream cur(fname.str(), ios::in);

				t1.setOrder(j);
				int str_amount = 0;
				while (!cur.eof())
				{
					str_amount++;
					if (str_amount % 200000 == 0)
						cout << "Or:" << j << " Sub:" << k << "\n";
					if (str_amount % 5000 == 0)
						cout << str_amount << " ";
					getline(cur, s);
					if (s.length() > 0)
					{
						istringstream iss;
						iss.str(s);
						iss >> tmp_s >> temp_val;


						t1.decompose(tmp_s, temp_val);
					
						ge.addTermRotation(j, t1);
					}
				}
				cur.close();
			}


		
			ge.printTermRotation(out_rot, j);
			if (j != max_order)
			out_rot << ",";
		

			std::cout << "Done!\n";
		}



		out_rot << "}";
		out_rot.close();

	}

	//mode 3.0
	void eval_excitations_numerical(int i,int antiferro_mode)
	{
		std::ofstream out_rot_excitation_numericals_single;
		AopXZRotateStorage rotate_excitation_storage;
		term t1;

		string s, tmp_s;
		double temp_val;

		if (antiferro_mode == 1)
			tmp_s = "_" + add_type;
		else
			tmp_s = "";
		fname.str("");
		fname.precision(2);
		fname << fixed;
		fname << out_res << delim << (*points)[i] << delim << "excitations_rot_" << (*points)[i] << "_" << max_order << "_" << angle_start2 << "_" << angle_start << tmp_s << "_single.txt";
		out_rot_excitation_numericals_single.open(fname.str(), ios::out);
		out_rot_excitation_numericals_single.precision(10);
		out_rot_excitation_numericals_single << fixed;
		out_rot_excitation_numericals_single << "{";
		//init of array for numerical excitations files
		double angle_cur = angle_start, angle_cur2 = angle_start2;

		for (int j = min_order; j <= max_order; j++)
		{
			rotate_excitation_storage.clearTerms();
			rotate_excitation_storage.set(m, size, false, angle_start, angle_finish, angle_step, angle_start2, angle_finish2, angle_step2, true);
			for (int k = min_op_amount; k <= j; k++)
			{
				cout << "SubOrder: " << k << "\n";
				fname.str("");
				fname << "d:\\Andrew\\Practice\\!!!_Last Set\\6.2 Collect Terms\\6.2 Collect Terms" << delim << inp_res << delim << (*points)[i] << "\\" << j << "_results_" << (*points)[i] << "_" << k << out_file_end;

				ifstream cur(fname.str(), ios::in);

				t1.setOrder(j);

				int str_amount = 0;
				while (!cur.eof())
				{
					str_amount++;
					if (str_amount % 200000 == 0)
						cout << "Or:" << j << " Sub:" << k << "\n";
					if (str_amount % 5000 == 0)
						cout << str_amount << " ";
					getline(cur, s);
					if (s.length() > 0)
					{
						istringstream iss;
						iss.str(s);
						iss >> tmp_s >> temp_val;


						t1.decompose(tmp_s, temp_val);
						
							if (t1.len > 0)//for rotate excitations
								if (antiferro_mode == 1)
								{
									rotate_excitation_storage.ConvertTo2TypesNumerical(t1, 0);
									rotate_excitation_storage.ConvertTo2TypesNumerical(t1, 1);
								}
								else
									rotate_excitation_storage.ConvertToBilinear(t1);
							break;
					}
				}
				cur.close();
			}
			rotate_excitation_storage.print_numerical_single_term_2_types(out_rot_excitation_numericals_single);
			if (j != max_order) {
				out_rot_excitation_numericals_single << ",";
			}
		}
		out_rot_excitation_numericals_single << "}";
		out_rot_excitation_numericals_single.close();
	}

	//mode==4
	void eval_energy_rotate_antiferro(int i)
	{
		groundEnergy ge;
		term t1;
		string s, tmp_s;
		double temp_val;

		ge.set(-0.5, true);

		std::ofstream out_energy_rot_antiferro;
		fname.str("");
		fname << out_res << delim << "energy_rot_" << add_type << "_" << (*points)[i] << "_" << max_order << ".txt";
		out_energy_rot_antiferro.open(fname.str(), ios::out);
		out_energy_rot_antiferro.precision(10);
		out_energy_rot_antiferro << fixed;
		out_energy_rot_antiferro << "{";

		for (int j = min_order; j <= max_order; j++)
		{
			ge.clearTerms();
			for (int k = min_op_amount; k <= j; k++)
			{
				cout << "SubOrder: " << k << "\n";
				fname.str("");
				fname << "d:\\Andrew\\Practice\\!!!_Last Set\\6.2 Collect Terms\\6.2 Collect Terms" << delim << inp_res << delim << (*points)[i] << "\\" << j << "_results_" << (*points)[i] << "_" << k << out_file_end;
				ifstream cur(fname.str(), ios::in);

				t1.setOrder(j);
				int str_amount = 0;
				while (!cur.eof())
				{
					str_amount++;
					if (str_amount % 200000 == 0)
						cout << "Or:" << j << " Sub:" << k << "\n";
					if (str_amount % 5000 == 0)
						cout << str_amount << " ";
					getline(cur, s);
					if (s.length() > 0)
					{
						istringstream iss;
						iss.str(s);
						iss >> tmp_s >> temp_val;

						t1.decompose(tmp_s, temp_val);

						ge.addTermRotationZeroPi2Angles(j, t1);
					}
				}
				cur.close();
			}
			ge.printTermRotation(out_energy_rot_antiferro, j, true);
			if (j != max_order)
				out_energy_rot_antiferro << ",";

		}
		out_energy_rot_antiferro << "}";
		out_energy_rot_antiferro.close();
	}




};

int main(int argc, char * argv[])
{
	cout << min_op_amount << "\n";
	vector<string> points;
	std::ostringstream fname;
	fname << config_dir << delim << "config.txt";
	ifstream config(fname.str(), ios::in);
	int num_points, min_order, max_order;
	string tmp;
	int a_amount,mode,analytical_mode,antiferro_mode;
	double angle_start, angle_step, angle_finish;
	double angle_start2, angle_step2, angle_finish2;

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
	config >> antiferro_mode;//0-false; 1 - true;

	//if there are no command line arguments read angles from file
	if (argc < 2)
	{
		getline(config, tmp);
		getline(config, tmp);
		config >> angle_start >> angle_step >> angle_finish;
		getline(config, tmp);
		getline(config, tmp);
		config >> angle_start2 >> angle_step2 >> angle_finish2;
	}
	else
	{
		cout << argv[1] << " " << argv[2]<<endl;
		angle_start = atof(argv[1]);
		angle_start2 = atof(argv[2]);
		angle_step = 100;
		angle_step2 = 100;
		angle_finish = angle_start - 1;
		angle_finish2 = angle_start2 - 1;
	}
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

	eval_results my_res(points,m,size);
	my_res.angle_start = angle_start;
	my_res.angle_start2 = angle_start2;
	my_res.angle_finish = angle_finish;
	my_res.angle_finish2 = angle_finish2;
	my_res.angle_step = angle_step;
	my_res.angle_step2 = angle_step2;

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
		
		
		//mode:
		//1 - energy + rotation
		//3 - excitations + rotation
		//4 - energy antiferro + rotation
		
		switch (mode)
		{
			case 1:
				my_res.eval_energy_rotate(i);
				break;
			case 3:
				my_res.eval_excitations_numerical(i, 1);
				break;
			case 4:
				my_res.eval_energy_rotate_antiferro(i);
				break;
		}
	}

	return 0;
}

