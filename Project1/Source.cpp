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

int Gauss(vector<vector<double>>& A, vector<double>& b, int N)
{
	// Гаусс
   // приведение к треугольному виду
	double t;
	for (int k = 0; k < N - 1; k++)
	{
		for (int i = k + 1; i < N; i++)
		{
			t = A[i][k] / A[k][k];
			b[i] -= t * b[k];
			for (int j = k + 1; j < N; j++)
			{
				A[i][j] -= t * A[k][j];
			}
		}
	}
	b[N - 1] /= A[N - 1][N - 1];
	// решение СЛАУ с треугольной матрицей
	for (int k = N - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < N; j++)
		{
			sum += A[k][j] * b[j];
		}
		b[k] = (b[k] - sum) / A[k][k];
	}

	return 0;
}


void main()
{
	//-----------------------------Init------------------------------------------
	ofstream ofp;
	bool write_all_results = 1;
 	setlocale(0, ""); // установка русского языка для вывода
	FEM_electro task;
	vector<double> V, V_tmp;
	string path = "test_for_me_2\\";
	
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
	int true_vel_id = 0;
	double true_I = abs(task.VELs[0].J[0]);
	//task.DirectTask();
	task.DirectTask(0, 1.0, true);
	task.V_in_rec(V, true_vel_id);
	for (int i = 0; i < V.size(); i++)
	{
		task.receivers[i].V_true = true_I * V[i];
		//task.receivers[i].V_true = 1.01 * true_I *  V[i]; // добавить шум 1% во все приемники
	}
	//task.receivers[4].V_true *= 1.1; // добавить шум 10%

	ifstream in(path + "Receiver.txt");
	ofp.open(path + "Receiver_V.txt");
	ofp << task.receivers.size() << endl;
	for (int i = 0; i < task.receivers.size(); i++)
	{
		/*ofp << task.receivers[i].x[0] << " " << task.receivers[i].y[0] << " "
			<< task.receivers[i].x[1] << " " << task.receivers[i].y[1] << " "
			<< task.receivers[i].V_true << endl;*/
		ofp << task.GetR(i, true_vel_id,0) << '\t' << task.GetR(i, true_vel_id, 1) << '\t' << task.receivers[i].V_true << endl;
	}
	ofp.close();
	//--------------------------------------------------------------------------

	/*for (int i = 0; i < V.size(); i++)
	{
		cout << " r " << task.receivers[i].r << '\t' << "num " << V[i] << '\t' << " analitic " << task.VELs[0].J[0] / task.mesh.sigma[0] / (sqrt(task.receivers[i].r * task.receivers[i].r + task.receivers[i].z * task.receivers[i].z)) << endl;
	}*/

	// -------------------------обратная задача---------------------------------
	// ---------------------------поиск сигмы-----------------------------------
	/*ofp.open(path + "history_func.txt");
	ofp << "iter " << '\t' << "sigma" << '\t' << "func" << endl;
	// начальные значения сигм
	double sigma_cur = 0.1; // 0.1
	int sigma_id = 0;
	double h;
	double eps = 1e-14;
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
	double prev_func = 1000;
	// функционал: SUM (1/V_true(V-Vtrue))^2 по всем приемникам
	for (int iter = 0; iter < max_iter && func_now > eps && abs(prev_func - func_now)/ func_now > 1e-8; iter++)
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
		
			// без весов (!в task.GetInverseFunc() тоже менять!)
			//a += ((V[i] - task.receivers[i].V) / h / 1.0) * ((V[i] - task.receivers[i].V) / h / 1.0);
			//b += -(task.receivers[i].V - task.receivers[i].V_true) * (V[i] - task.receivers[i].V) / h / 1.0;

		}
		
		// сигма_след = сигма + b1 / а11
		sigma_cur = sigma_cur + b / a;

		task.DirectTask(sigma_id, sigma_cur);
		task.V_in_rec(V);
		for (int i = 0; i < V.size(); i++)
		{
			task.receivers[i].V = V[i];
		}
		prev_func = func_now;
		func_now = task.GetInverseFunc();
		cout << "iter " << iter + 1 << " sigma " << sigma_cur << " func " << func_now << " d_func " << abs(prev_func - func_now) / func_now << endl;
		ofp << iter + 1 << '\t' << sigma_cur << '\t' << func_now << '\t' << abs(prev_func - func_now) / func_now << endl;
	}
	ofp.close();*/

	// ---------------------------поиск силы тока-----------------------------------
	/*ofp.open(path + "history_func_I.txt");
	ofp << "iter " << '\t' << "I" << '\t' << "func" << endl;
	// начальные значения
	double I_cur = 6;
	double h;
	double eps = 1e-14;
	int max_iter = 30;
	double func_now = 100;
	double a, b;
	task.DirectTask(0, I_cur, true);
	task.V_in_rec(V);
	for (int i = 0; i < V.size(); i++)
	{
		task.receivers[i].V = V[i];
	}
	func_now = task.GetInverseFunc();
	cout << "iter " << 0 << " I " << I_cur << " func " << func_now << endl;
	ofp << 0 << '\t' << I_cur << '\t' << func_now << endl;
	double prev_func = 1000;
	// функционал: SUM (1/V_true(V-Vtrue))^2 по всем приемникам
	for (int iter = 0; iter < max_iter && func_now > eps && abs(prev_func - func_now) / func_now > 1e-8; iter++)
	{
		h = 0.05 * I_cur;
		task.DirectTask(0, I_cur + h, true);
		task.V_in_rec(V);

		// ищем 1 сигму, поэтому СЛАУ а11 шаг =b1, где а11 = SUM ( (1/V_true(V(sigma + h) - V(sigma))/h)^2), b1 = SUM ( (1/V_true)^2*(V(sigma + h) - V(sigma))/h*(V-Vtrue))
		a = 0;
		b = 0;
		for (int i = 0; i < V.size(); i++)
		{
			a += ((V[i] - task.receivers[i].V) / h / task.receivers[i].V_true) * ((V[i] - task.receivers[i].V) / h / task.receivers[i].V_true);
			b += -(task.receivers[i].V - task.receivers[i].V_true) * (V[i] - task.receivers[i].V) / h / task.receivers[i].V_true / task.receivers[i].V_true;

			// без весов (!в task.GetInverseFunc() тоже менять!)
			//a += ((V[i] - task.receivers[i].V) / h / 1.0) * ((V[i] - task.receivers[i].V) / h / 1.0);
			//b += -(task.receivers[i].V - task.receivers[i].V_true) * (V[i] - task.receivers[i].V) / h / 1.0;

		}

		// сигма_след = сигма + b1 / а11
		I_cur = I_cur + b / a;

		task.DirectTask(0, I_cur, true);
		task.V_in_rec(V);
		for (int i = 0; i < V.size(); i++)
		{
			task.receivers[i].V = V[i];
		}
		prev_func = func_now;
		func_now = task.GetInverseFunc();
		cout << "iter " << iter + 1 << " I " << I_cur << " func " << func_now << " d_func " << abs(prev_func - func_now) / func_now << endl;
		ofp << iter + 1 << '\t' << I_cur << '\t' << func_now << '\t' << abs(prev_func - func_now) / func_now << endl;
	}*/

	// ---------------------------поиск местоположения-----------------------------------
	ofp.open(path + "history_func_I_2.txt");
	ofp << "iter " << '\t' << "alpha_loc" << '\t';
	for (int i = 0; i < 10; i++)
	{
		ofp << "I" << i << '\t';
	}
	ofp << "func" << '\t' << "d_func" << endl;
	// начальные значения
	vector<double> I_cur =
	{
	9.5238,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0
	};
	
	double h;
	double eps = 1e-14;
	int max_iter = 5;
	double func_now = 100, func_save, func_alpha;
	double alpha = 0, alpha_loc;

	vector<vector<double>> dV; // dV[i] - V с приращением по iому параметру
	dV.resize(10);
	for (int i = 0; i < 10; i++)
	{
		dV[i].resize(V.size());
	}
	vector<vector<double>> V_one_source; // V_one_source[i] - V от iого источника
	V_one_source.resize(10);
	for (int i = 0; i < 10; i++)
	{
		V_one_source[i].resize(V.size());
	}
	vector<vector<double>> A_orig;
	A_orig.resize(10);
	for (int i = 0; i < 10; i++)
	{
		A_orig[i].resize(10);
	}
	vector<double> B_orig(10);
	vector<vector<double>> A;
	A.resize(10);
	for (int i = 0; i < 10; i++)
	{
		A[i].resize(10);
	}
	vector<double> B(10), B_save(10), V_save;

	fill(V.begin(), V.end(), 0);
	for (int i_vel = 0; i_vel < 10; i_vel++)
	{
		task.DirectTask(0, I_cur[i_vel], true);
		task.V_in_rec(V_tmp, i_vel);
		for (int i = 0; i < V.size(); i++)
		{
			V[i] += V_tmp[i];
			V_one_source[i_vel][i] = V_tmp[i];
		}
	}
	for (int i = 0; i < V.size(); i++)
	{
		task.receivers[i].V = V[i];
	}
	
	//func_now = task.GetInverseFunc(V, 0, I_cur);
	func_now = task.GetInverseFunc();

	cout << "iter " << 0 << " func " << func_now << endl;
	for (int i = 0; i < 10; i++)
	{
		cout << I_cur[i] << endl;
	}
	ofp << 0 << '\t' << "   " << '\t';
	for (int i = 0; i < 10; i++)
	{
		ofp << I_cur[i] << '\t';
	} 
	ofp << func_now << endl;
	double prev_func = 1000, prev_func_alpha;
	// функционал: SUM (1/V_true(V-Vtrue))^2 по всем приемникам
	for (int iter = 0; iter < max_iter && func_now > eps && abs(prev_func - func_now) / func_now > 1e-8; iter++)
	{
		// запомним значения с приращением
		for (int i_param = 0; i_param < 10; i_param++)
		{
			fill(dV[i_param].begin(), dV[i_param].end(), 0);
			for (int i_vel = 0; i_vel < 10; i_vel++)
			{
				if (i_param == i_vel)
				{
					h = 0.05 * I_cur[i_vel];
					task.DirectTask(0, I_cur[i_vel] + h, true);
					task.V_in_rec(V_tmp, i_vel);
					for (int i = 0; i < V.size(); i++)
					{
						dV[i_param][i] += V_tmp[i];
					}
				}
				else
				{
					h = 0;
					for (int i = 0; i < V.size(); i++)
					{
						dV[i_param][i] += V_one_source[i_vel][i];
					}
				}
				
			}
		}
		// значение без приращения есть в V
		
		// сборка матрицы
		for (int i = 0; i < 10; i++)
		{
			B_orig[i] = 0;
			for (int j = 0; j < 10; j++)
			{
				A_orig[i][j] = 0;
				for (int k = 0; k < V.size(); k++)
				{
					A_orig[i][j] +=
						(1.0 / (task.receivers[k].V_true * task.receivers[k].V_true))  // w_k
						* ((dV[i][k] - V[k]) / (0.05 * I_cur[i]))  // dVk/dI_i
						* ((dV[j][k] - V[k]) / (0.05 * I_cur[j])); // dVk/dI_j
					B_orig[i] +=
						(-1.0 / (task.receivers[k].V_true * task.receivers[k].V_true))  // w_k
						* ((dV[i][k] - V[k]) / (0.05 * I_cur[i]))  // dVk/dI_i
						* (V[k] - task.receivers[k].V_true); // (V(I) - V*)
				}
			}
		}
		
		// регуляризация
		alpha_loc = alpha;
		A = A_orig;
		B = B_orig;
		for (int i = 0; i < 10; i++)
		{
			A[i][i] += alpha_loc;
		}
		Gauss(A, B, 10);
		// подстановка новой I
		fill(V.begin(), V.end(), 0);
		for (int i_vel = 0; i_vel < 10; i_vel++)
		{
			task.DirectTask(0, I_cur[i_vel] + B[i_vel], true);
			task.V_in_rec(V_tmp, i_vel);
			for (int i = 0; i < V.size(); i++)
			{
				V[i] += V_tmp[i];
				V_one_source[i_vel][i] = V_tmp[i];

			}
		}
		func_alpha = task.GetInverseFunc(V, alpha_loc, I_cur, B);
		prev_func_alpha = 100;
		// подбор альфы
		do
		{
			ofp << "   " << '\t' << alpha_loc << '\t';
			for (int i = 0; i < 10; i++)
			{
				ofp << I_cur[i] + B[i] << '\t';
			}
			ofp << func_alpha << '\t' << abs(prev_func_alpha - func_alpha) / func_alpha << endl;
			cout << "   " << '\t' << alpha_loc << '\t';
			for (int i = 0; i < 10; i++)
			{
				cout << I_cur[i] + B[i] << '\t';
			}
			cout << func_alpha << '\t' << abs(prev_func_alpha - func_alpha) / func_alpha << endl;
			B_save = B;
			V_save = V;
			func_save = func_alpha;
			alpha_loc *= 10;
			A = A_orig;
			B = B_orig;
			for (int i = 0; i < 10; i++)
			{
				A[i][i] += alpha_loc;
			}
			Gauss(A, B, 10);
			// подстановка новой I
			fill(V.begin(), V.end(), 0);
			for (int i_vel = 0; i_vel < 10; i_vel++)
			{
				task.DirectTask(0, I_cur[i_vel] + B[i_vel], true);
				task.V_in_rec(V_tmp, i_vel);
				for (int i = 0; i < V.size(); i++)
				{
					V[i] += V_tmp[i];
					V_one_source[i_vel][i] = V_tmp[i];

				}
			}
			prev_func_alpha = func_alpha;
			func_alpha = task.GetInverseFunc(V, alpha_loc, I_cur, B);
		} while (prev_func_alpha > func_alpha && alpha_loc < 1e-5);
		ofp << "   " << '\t' << alpha_loc << '\t';
		for (int i = 0; i < 10; i++)
		{
			ofp << I_cur[i] + B[i] << '\t';
		}
		ofp << func_alpha << '\t' << abs(prev_func_alpha - func_alpha) / func_alpha << endl;
		cout << "   " << '\t' << alpha_loc << '\t';
		for (int i = 0; i < 10; i++)
		{
			cout << I_cur[i] + B[i] << '\t';
		}
		cout << func_alpha << '\t' << abs(prev_func_alpha - func_alpha) / func_alpha << endl;

		// сдвиг сил тока
		for (int i_vel = 0; i_vel < 10; i_vel++)
		{
			I_cur[i_vel] += B_save[i_vel];
		}
		for (int i = 0; i < V.size(); i++)
		{
			task.receivers[i].V = V_save[i];
		}
		
		prev_func = func_now;
		func_now = func_save;

		cout << "iter " << iter + 1;
		for (int i = 0; i < 10; i++)
		{
			cout << I_cur[i] << '\t';
		}
		cout << " func " << func_now << " d_func " << abs(prev_func - func_now) / func_now << endl;
		ofp << iter + 1 << '\t' << "   " << '\t';
		for (int i = 0; i < 10; i++)
		{
			ofp << I_cur[i] << '\t';
		} 
		ofp << func_now << '\t' << abs(prev_func - func_now) / func_now << endl;
	}

	return;
}
