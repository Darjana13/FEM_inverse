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

	//-----------------------------Solve----------------------------------------
	task.InitTask(path);
	task.DirectTask();
	task.V_in_rec(V);

	task.PrintField("fieldV.txt");

	return;
}
