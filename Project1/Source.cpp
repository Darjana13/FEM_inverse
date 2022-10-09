#include "includes.h"
#include "FEM.h"

int GetPoints(vector<pair<double, double>> &point, string path)
{
	int i;

	ifstream inf(path + "Point");
	if (!inf)
	{
		cout << "No file " << path << "Point" << endl;
		return 1;
	}
	int kol;
	inf >> kol;
	point.resize(kol);
	double x, y;
	for (i = 0; i < kol; i++)
	{
		inf >> x >> y;
		point[i] = make_pair(x, y);
	}

	inf.close();
	inf.clear();

	return 0;
}

void main()
{
	//-----------------------------Init------------------------------------------
	ofstream ofp;
	bool write_all_results = 1;
 	setlocale(0, ""); // установка русского языка для вывода
	FEM_electro task;
	vector<double> V;
	string path = "test\\";
	//vector<pair<double, double>> points;
	//GetPoints(points, path);
	//----------------------------------------------------------------------------

	/*task.InitTask(path);
	vector<double> q_ist(task.mesh.kol_nodes);
	ofp.open("true_field.txt");
	for (int i = 0; i < task.mesh.kol_nodes; i++)
	{
		q_ist[i] = 1. / 0.1 / (sqrt(task.mesh.nodes[i].r * task.mesh.nodes[i].r + task.mesh.nodes[i].z * task.mesh.nodes[i].z));
		if (q_ist[i] != q_ist[i])
			q_ist[i] = 1. / 0.1 / 0.0014;

		ofp << task.mesh.nodes[i].r << " " << task.mesh.nodes[i].z << " " << q_ist[i] << endl;
	}*/

	//-----------------------------Solve----------------------------------------
	task.InitTask(path);
	task.DirectTask();
	task.V_in_rec(V);

	task.PrintField("fieldV.txt");
	for (int i = 0; i < V.size(); i++)
	{
		cout << " r " << task.receivers[i].r << '\t' << "num " << V[i] << '\t' << " analitic " << task.VELs[0].J[0] / task.mesh.sigma[0] / (sqrt(task.receivers[i].r * task.receivers[i].r + task.receivers[i].z * task.receivers[i].z)) << endl;
	}
	return;
}
