#pragma once
#include "stdafx.h"
#ifndef __PSEUDOSPIN_TYPES_H__
#define __PSEUDOSPIN_TYPES_H__

const int Types = 3; //amount of different types of rotation angels


class pseudospin_types
{
public:
	static int getType(int n, int zeroType);
	static std::string getAngleName(int type);
	static void fillMatrix();
	static void setOrder(int NNN);
private:
	static const std::string angleNames[Types]; //names of different types of rotation angels
	static int **matrix;
	static int NNN;
	static int size();
};



#endif  //__PSEUDOSPIN_TYPES_H__