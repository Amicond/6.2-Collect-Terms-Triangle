#include "stdafx.h"
#include "a_op_xz_rotate.h"



//void AopXZRotate::add_analytic_operator(term &t, int type, int n){//type: 0 - S_plus, 1 - S_minus
//	aop.node.push_back(t.nums[n]); //set first node-number
//
//	if (type == 0){ //select Sp' -term for first a-operator
//		aop.names.push_back('p');
//		switch (t.ops[n])	{
//		case 'z':
//			aop.coeff *= -0.5;
//			aop.coeff *= 2;
//			trc.incCoeff(2, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
//			trc.incCoeff(3, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
//			break;
//		case 'p':
//			//current.aop.coeff *= 1;
//			trc.incCoeff(2, 2);//Cos^2[beta/2]
//			break;
//		case 'm':
//			aop.coeff *= -1;
//			trc.incCoeff(3, 2);//Sin^2[beta/2]
//			break;
//		}
//	}
//	else{ //select Sm-term for first a-operator
//		aop.names.push_back('m');
//		switch (t.ops[n])	{
//		case 'z':
//			aop.coeff *= -0.5;
//			aop.coeff *= 2;
//			trc.incCoeff(2, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
//			trc.incCoeff(3, 1);//Sin[beta]=2*Cos[beta/2]*Sin[beta/2]
//			break;
//		case 'p':
//			aop.coeff *= -1;
//			trc.incCoeff(3, 2);//Sin^2[beta/2]
//			break;
//		case 'm':
//			//current.aop.coeff *= 1;
//			trc.incCoeff(2, 2);//Cos^2[beta/2]
//			break;
//		}
//	}
//}

bool AopXZRotate::operator==(const AopXZRotate TxzR2)const
{
	if (!(aop_c == TxzR2.aop_c))
		return false;
	else if (!(trc == TxzR2.trc))
		return false;
	return true;
}

bool AopXZRotate::operator<(const AopXZRotate TxzR2)const
{
	if (aop_c < TxzR2.aop_c) return true;
	else
		if (TxzR2.aop_c < aop_c)return false;
		else
			if (trc < TxzR2.trc) return true;
			else return false;
}

void AopXZRotate::add_trig(int type, int power, bool second ){
	if (second)
		trc2.incCoeff(type, power);
	else
		trc.incCoeff(type, power);

}

void AopXZRotate::clear()
{
	aop_c.clear();
	trc.clear();
	trc2.clear();
	sz_power = 0;
}
