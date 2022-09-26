#include "Table.h"

Table::Table(string path)
{
	ReadTable(path);
}

int Table::ReadTable(string path)
{
	ifstream inf;
	inf.open(path);
	inf >> N;
	x_all.resize(N);
	f_all.resize(N);
	q.resize(N);
	A.resize(N, vector<double>(N, 0.0));
	b.resize(N, 0.0);
	w.resize(N, 1.0);

	for (int i = 0; i < N; i++)
	{
		inf >> f_all[i] >> x_all[i];
	}
	return 0;
}

double Table::GetF_linear_interpolation(double _x)
{
	// find interval
	if (N == 0)
		return 0;
	if (_x < x_all[0])
		return f_all[0];
	if (_x > x_all[N - 1])
		return _x / (x_all[N - 1] * (1 - f_all[N - 1]) / f_all[N - 1] + _x);
	for (int i = 0; i < N - 1; i++)
	{
		if (x_all[i] <= _x && _x <= x_all[i + 1])
		{
			return f_all[i] + (f_all[i+1] - f_all[i]) / (x_all[i+1] - x_all[i]) * (_x - x_all[i]);
		}
	}
}

// do not work
/*double Table::GetF(double _x) // return mu_otnositel'noe
{
	double y = 0;
	double ksi = 0, h = 0;
	// find interval
	if (_x < x_all[0])
		return 0;
	if (_x > x_all[N-1])
		return _x / (x_all[N-1] * (1 - f_all[N-1])/ f_all[N - 1] + _x);
	for (int i = 0; i < N - 1; i++)
	{
		if (x_all[i] <= _x && _x <= x_all[i + 1])
		{
			ksi = getKsi(_x, x_all[i], x_all[i + 1]);
			h = x_all[i], x_all[i + 1];
			for (int j = 0; j < 4; j++)
			{		
				y += ((j % 4) % 2 ? h : 1) * ermit_func[j % 4](ksi);
			}
		}
	}


	return y;
}
*/
// do not work
/*int Table::Calc_q()
{
	double ksi = 0, h = 0;
	for (int i = 0; i < N; i++)
	{
		for (int k = 1; k < N - 1; k++)
		{
			ksi = getKsi(x_all[k], x_all[k], x_all[k + 1]);
			h = x_all[k], x_all[k + 1];
			for (int j = 0; j < N; j++)
			{// ????????		
				A[i][j] += w[j] *
					((i % 4) % 2 ? h : 1) * ermit_func[i % 4](ksi) *
					((j % 4) % 2 ? h : 1) * ermit_func[j % 4](ksi);
			}
			b[i] += w[i] * f_all[i] * ((i % 4) % 2 ? h : 1) * ermit_func[i % 4](ksi);
		}
	}
	Solve_Gauss(A, b, q);

	return 0;
}
*/