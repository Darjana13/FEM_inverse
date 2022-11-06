#include "FEM.h"

double FEM_electro::GetInverseFunc()
{
	double res = 0;
	for (int i = 0; i < receivers.size(); i++)
	{
		res += (receivers[i].V - receivers[i].V_true)* (receivers[i].V - receivers[i].V_true) / (receivers[i].V_true * receivers[i].V_true);
		
		// без весов (!в main тоже менять!)
		//res += (receivers[i].V - receivers[i].V_true) * (receivers[i].V - receivers[i].V_true);

		//cout << " V_true " << receivers[i].V_true << " V " << receivers[i].V << " res " << res << endl;
	}
	//cout << " res " << res << endl;
	return res;
}

double FEM_electro::GetInverseFunc(vector<double>& V, double alpha, vector<double>& I)
{
	double res = 0;
	for (int i = 0; i < receivers.size(); i++)
	{
		res += (V[i] - receivers[i].V_true) * (V[i] - receivers[i].V_true) / (receivers[i].V_true * receivers[i].V_true);
	}
	//cout << " res1 " << res << endl;

	if (alpha != 0)
	{
		for (int i = 0; i < I.size(); i++)
		{
			res += alpha * (I[i]*I[i]);
		}
	}
	//cout << " res2 " << res << endl;
	return res;
}

double FEM_electro::GetInverseFunc(vector<double>& V, double alpha, vector<double>& I, vector<double>& Ih)
{
	double res = 0;
	for (int i = 0; i < receivers.size(); i++)
	{
		res += (V[i] - receivers[i].V_true) * (V[i] - receivers[i].V_true) / (receivers[i].V_true * receivers[i].V_true);
	}
	//cout << " res1 " << res << endl;

	if (alpha != 0)
	{
		for (int i = 0; i < I.size(); i++)
		{
			res += alpha * ((I[i] + Ih[i]) * (I[i] + Ih[i]));
		}
	}
	//cout << " res2 " << res << endl;
	return res;
}
