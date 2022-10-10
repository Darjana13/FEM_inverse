#include "FEM.h"

double FEM_electro::GetInverseFunc()
{
	double res = 0;
	for (int i = 0; i < receivers.size(); i++)
	{
		res += (receivers[i].V - receivers[i].V_true)* (receivers[i].V - receivers[i].V_true) / receivers[i].V_true;

		cout << " V_true " << receivers[i].V_true << " V " << receivers[i].V << " res " << res << endl;
	}
	cout << " res " << res << endl;
	return res;
}