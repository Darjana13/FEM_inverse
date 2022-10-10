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
	string path = "test 3layers\\";
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
	task.InitTask(path); // создаем сетку и читаем данные

	//-------------------получаем синтетические данные--------------------------
	task.DirectTask();
	task.V_in_rec(V);
	for (int i = 0; i < V.size(); i++)
	{
		task.receivers[i].V_true = V[i];
	}
	//--------------------------------------------------------------------------

	/*for (int i = 0; i < V.size(); i++)
	{
		cout << " r " << task.receivers[i].r << '\t' << "num " << V[i] << '\t' << " analitic " << task.VELs[0].J[0] / task.mesh.sigma[0] / (sqrt(task.receivers[i].r * task.receivers[i].r + task.receivers[i].z * task.receivers[i].z)) << endl;
	}*/

	// -------------------------обратная задача---------------------------------

	ofp.open(path + "history_func.txt");
	ofp << "iter " << '\t' << "sigma" << '\t' << "func" << endl;
	// начальные значения сигм
	double sigma_cur = 0.01; // 0.1
	int sigma_id = 1;
	double h;
	double eps = 1e-10;
	int max_iter = 30;
	double func_now = 100;
	double a, b;
	task.DirectTask(0, sigma_cur);
	task.V_in_rec(V);
	for (int i = 0; i < V.size(); i++)
	{
		task.receivers[i].V = V[i];
	}
	func_now = task.GetInverseFunc();
	cout << "iter " << 0 << " sigma " << sigma_cur << " func " << func_now << endl;
	ofp << 0 << '\t' << sigma_cur << '\t' << func_now << endl;

	// функционал: SUM (1/V_true(V-Vtrue)^2) по всем приемникам
	for (int iter = 0; iter < max_iter && func_now > eps; iter++)
	{
		h = 0.05 * sigma_cur;
		task.DirectTask(sigma_id, sigma_cur + h);
		task.V_in_rec(V);

		// ищем 1 сигму, поэтому СЛАУ а11 шаг =b1, где а11 = SUM ( (1/V_true(V(sigma + h) - V(sigma))/h)^2), b1 = SUM ( (1/V_true)^2*(V(sigma + h) - V(sigma))/h*(V-Vtrue))
		a = 0;
		b = 0;
		for (int i = 0; i < V.size(); i++)
		{
			a += ( (V[i] - task.receivers[i].V) / h / task.receivers[i].V_true)* ((V[i] - task.receivers[i].V) / h / task.receivers[i].V_true);
			b += -(task.receivers[i].V - task.receivers[i].V_true) * (V[i] - task.receivers[i].V) / h / task.receivers[i].V_true / task.receivers[i].V_true;
		}
		
		// сигма_след = сигма + b1 / а11
		sigma_cur = sigma_cur + b / a;

		task.DirectTask(sigma_id, sigma_cur);
		task.V_in_rec(V);
		for (int i = 0; i < V.size(); i++)
		{
			task.receivers[i].V = V[i];
		}
		func_now = task.GetInverseFunc();
		cout << "iter " << iter + 1 << " sigma " << sigma_cur << " func " << func_now << endl;
		ofp << iter + 1 << '\t' << sigma_cur << '\t' << func_now << endl;
	}
	
	return;
}
