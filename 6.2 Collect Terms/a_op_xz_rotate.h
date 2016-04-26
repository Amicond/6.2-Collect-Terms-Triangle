#include "stdafx.h"
#ifndef __AOP_XZ_ROTATE_H__
#define  __AOP_XZ_ROTATE_H__
#include "term.h"
#include "a_op_couple.h"
#include "trig_coefficients.h"


class AopXZRotate {
	//bilinear a-operator Term + Trig Function according to spin rotation
public:
	a_op_couple aop_c;
	trig_coefficients trc,trc2;
	int sz_power; //power of sz-factor (sz is equal to 1/2 or -1/2)
//Methods
	//void add_analytic_operator(term &t, int type, int n);//type: 0 - S_plus, 1 - S_minus
	AopXZRotate()
	{

	}

	bool operator==(const AopXZRotate TxzR2)const;
	
	bool operator<(const AopXZRotate TxzR2)const;

	void add_trig(int type, int power, bool second = false);


	/*void setHash();

	size_t getHash() const;*/
	
	void clear();
};

#endif //end of __AOP_XZ_ROTATE_H__