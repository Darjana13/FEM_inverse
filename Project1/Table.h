#pragma once
#include "includes.h"
#include "cgm.h"

class Table
{
	// для сплайна; не работает
	//int Calc_q();
	//double GetF(double _x);

public:
	vector<double> x_all, f_all;
	vector<vector<double>> A;
	vector<double> q, b;
	vector<double> w;
	int N;

	int ReadTable(string path);

	double GetF_linear_interpolation(double _x);

	Table(string path);
	Table(vector<double> x, vector<double> f)
	{
		x_all = x;
		f_all = f;
		N = x.size();
		q.resize(N);
		A.resize(N, vector<double>(N, 0.0));
		b.resize(N, 0.0);
		w.resize(N, 1.0);
	}
	Table()
	{
		N = 0;
	}

};