#include "stdafx.h"
#ifndef __A_OP_XZ_ROTATE_STORAGE__
#define __A_OP_XZ_ROTATE_STORAGE__
#include "a_op_xz_rotate.h"
#include "cos.h"

//storage class that contain all current elements

namespace std{
	template <>
	struct hash<a_op_couple>
	{
		std::size_t operator()(const a_op_couple& k) const;
	};
};

namespace std{
	template <>
	struct hash < AopXZRotate >
	{
		std::size_t operator()(const AopXZRotate& k) const;
	};
}

class AopXZRotateStorage{
public:
	//members
	std::unordered_map<AopXZRotate, double> storage; //for terms with trigonometry coeffs
	std::vector<std::unordered_map<a_op_couple, double>> storage_numerical;//for terms with numeric values only
	std::vector<std::vector<std::unordered_map<a_op_couple, double>>> storage_numerical_2_angle;//for terms with numeric values only
	int **matrix; //matrix vith node numbers
	int matrix_size;
	bool analytical_mode; //true- analytical, false -numerical;
	double angle_start, angle_start_2;
	double angle_finish, angle_finish_2;
	double angle_step, angle_step_2;
	double sz;
	static int arr[4][2];
	//methods
	//void set(int **Matrix, int Matrix_size, bool Mode, double Angle_start, double Angle_finish, double Angle_step, double Sz = -0.5);

	void set(int **Matrix, int Matrix_size, bool Mode, double Angle_start, double Angle_finish, double Angle_step, double Angle_start_2=0, double Angle_finish_2=-1, double Angle_step_2=1, double Sz = -0.5);

	void clearTerms();

	std::vector<std::pair<int, int>> generate_pairs(int n);

	void add_numerical(AopXZRotate &current);

	void add_numerical_2_angles(AopXZRotate &current);

	void add_operator(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos);

	void add_operator_pi_zero(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos);

	void add_operator_zero_pi_2_angles(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos);

	void add_operator_antiferro(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos);

	void add_operator_2_sublattices(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos);

	void add_operator_ZY2_sublattices(AopXZRotate &current, term &t, int num, int operator_type, int new_operator_pos);

	void ConvertToBilinear(term t);

	void ConvertToBilinearPiZero(term t);

	void ConvertToBilinearZeroPi2Angles(term t);

	void ConvertToBilinearAntiFerro(term t);

	void ConvertToBilinear2sublattices(term t);

	void ConvertToBilinearZY2sublattices(term t);

	void print(std::ofstream &F);

	void print_numerical(std::vector<std::ofstream*> &Fs);

	void print_numerical_2_angles(std::vector<std::vector<std::ofstream*>> &Fs);
};

#endif //__A_OP_XZ_ROTATE_STORAGE__