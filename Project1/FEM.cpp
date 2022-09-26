#include "FEM.h"

int FEM_electro::InitTask(string path)
{
	mesh.ReadMesh(path);
	ReadData(path);
	basis = { quadratic1, quadratic2, quadratic3 };
	slau.GeneratePortret(mesh);
	q.resize(mesh.kol_nodes);
	b_loc.resize(9);
	A_loc.resize(9);
	for (int i = 0; i < 9; i++)
	{
		A_loc[i].resize(9);
	}
	return 0;
}

int FEM_electro::ReadData(string path)
{
	ifstream in(path + "Receiver.txt");
	
	in >> kol_rec;
	receivers.resize(kol_rec);
	for (int i = 0; i < kol_rec; i++)
	{
		in >> receivers[i].r >> receivers[i].z >> receivers[i].V;
	}
	in.close();

	in.open(path + "VEL.txt");
	in >> kol_vel;
	VELs.resize(kol_vel);
	int cur_el = 0;
	double r1, r2, z1, z2;
	for (int i = 0; i < kol_vel; i++)
	{
		cur_el = -1;
		in >> VELs[i].r[0] >> VELs[i].z[0] >> VELs[i].r[1] >> VELs[i].z[1] >> VELs[i].J[0];
		VELs[i].J[1] = -VELs[i].J[0];
	}
	in.close();

	return 0;
}

int FEM_electro::GetG_Loc_biquadratic(int el_id) // получение локальной G
{
	double lambda = mesh.sigma[mesh.elements[el_id].material];
	double
		hr = mesh.nodes[mesh.elements[el_id].node_loc[2]].r - mesh.nodes[mesh.elements[el_id].node_loc[0]].r,
		hz = mesh.nodes[mesh.elements[el_id].node_loc[6]].z - mesh.nodes[mesh.elements[el_id].node_loc[0]].z;
	double r_cur = mesh.nodes[mesh.elements[el_id].node_loc[0]].r;

	int r_i, z_i, r_j, z_j;

	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			r_i = mu(i);
			r_j = mu(j);
			z_i = nu(i);
			z_j = nu(j);
			A_loc[i][j] = lambda * 
				((r_cur / 3./ hr * G_quadratic_r_1[r_i][r_j] + 1. / 6. * G_quadratic_r_2[r_i][r_j])*hz* M_quadratic[z_i][z_j]
				+ (hr*r_cur/30.*M_quadratic_r_1[r_i][r_j] + hr*hr/60.* M_quadratic_r_2[r_i][r_j])*1./hz*G_quadratic[z_i][z_j]);
		
			if (A_loc[i][j] != A_loc[i][j])
				cout << "nan" << endl;
		}
	}
	
	return 0;
}

double FEM_electro::GetUTest(double r, double z)
{
	return z*z*z+r*r*r;
}

double FEM_electro::GetF(double r, double z)
{
	//return -6*z-9*r;
	for (int i = 0; i < kol_vel; i++)
	{
		if (VELs[i].r[0] == r && VELs[i].z[0] == z)
		{
			cout << r << " " << z << " " << VELs[i].J[0] << endl;
			return VELs[i].J[0];
		}
		if (VELs[i].r[1] == r && VELs[i].z[1] == z)
			//return VELs[i].J[1];
			return 0;
	}
	return 0;
}

int FEM_electro::Getb_Loc_biquadratic(int el_id) // получение локального b
{
	vector<double> f(9);
	bool not_null = 0;
	for (int i = 0; i < 9; i++)
	{
		f[i] = GetF(mesh.nodes[mesh.elements[el_id].node_loc[i]].r, mesh.nodes[mesh.elements[el_id].node_loc[i]].z);
		if (f[i] != 0)
		{
			not_null = 1;
		}
	}

	if (not_null == 1)
	{
		double
			hr = mesh.nodes[mesh.elements[el_id].node_loc[2]].r - mesh.nodes[mesh.elements[el_id].node_loc[0]].r,
			hz = mesh.nodes[mesh.elements[el_id].node_loc[6]].z - mesh.nodes[mesh.elements[el_id].node_loc[0]].z;
		double r_cur = mesh.nodes[mesh.elements[el_id].node_loc[0]].r;

		int r_i, z_i, r_j, z_j;

		vector<vector<double>> M(9, vector<double>(9));
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				r_i = mu(i);
				z_i = nu(i);
				r_j = mu(j);
				z_j = nu(j);
				M[i][j] = (hr * r_cur / 30. * M_quadratic_r_1[r_i][r_j] + hr * hr / 60. * M_quadratic_r_2[r_i][r_j]) * hz * M_quadratic[z_i][z_j];
			}
		}


		for (int i = 0; i < 9; i++)
		{
			b_loc[i] = 0;
			r_i = mu(i);
			z_i = nu(i);
			for (int j = 0; j < 9; j++)
			{
				r_j = mu(j);
				z_j = nu(j);
				b_loc[i] += f[j] * (hr * r_cur / 30. * M_quadratic_r_1[r_i][r_j] + hr * hr / 60. * M_quadratic_r_2[r_i][r_j]) * hz * M_quadratic[z_i][z_j];
			}
		}
	}
	else
	{
		for (int i = 0; i < 9; i++)
		{
			b_loc[i] = 0;
		}
	}

	return 0;
}

int FEM_electro::DirectTask()
{
	slau.ClearValues();
	cout << "start make matrix " << endl;
	for (int i = 0; i < mesh.kol_elements; i++)
	{
		Getb_Loc_biquadratic(i);
		GetG_Loc_biquadratic(i);

		slau.AddLocal(A_loc, b_loc, i, mesh);
		if(i % 100 == 0)
			cout << "el " << i << " from " << mesh.kol_elements << endl;
	}
	//slau.SetS1(mesh);
	slau.SetS1_null(mesh);
	cout << "Set S1 " << endl;

	cgm_solver solver(slau.ia, slau.ja, slau.di, slau.au, slau.b);
	solver.llt_preconditioning(q);

	/*Solver solver(slau.ia, slau.ja, slau.di, slau.au, slau.al, slau.b);
	cout << "start solve " << endl;
	solver.CGM_LU();
	solver.getx0(q);*/

	/*Solver check(slau.ia, slau.ja, slau.di, slau.au, slau.al, slau.b);
	vector<double> b_ist(mesh.kol_nodes);
	solver.A.Ax(q, b_ist);
	for (int i = 0; i < mesh.kol_nodes; i++)
	{
		if (abs(slau.b[i] - b_ist[i]) > 1e-5)
			cout << " i " << i << " b[i] " << slau.b[i] << " Aq[i] " << b_ist[i] << endl;
	}*/


	return 0;
}

int FEM_electro::V_in_rec(vector<double> &V)
{
	/*double res = 0;
	for (int i = 0; i < mesh.kol_nodes; i++)
	{
		double f = GetUTest(mesh.nodes[i].r, mesh.nodes[i].z);
		if(f != 0)
			res += abs(q[i] - f) / abs(f);
		else
			res += abs(q[i] - f);
	}
	cout << "res " << res / mesh.kol_nodes << endl;*/

	V.resize(kol_rec);
	for(int i = 0; i < kol_rec; i++)
	{
		V[i] = V_in_point(receivers[i].r, receivers[i].z);
	}
	return 0;
}

int FEM_electro::PrintField(string filename)
{
	ofstream out(filename);
	for (int i = 0; i < mesh.kol_nodes; i+=1)
	{
		out << mesh.nodes[i].r << " " << mesh.nodes[i].z << " " << q[i] << endl;
		/*if (q[i] < 0)
		{
			cout << " q[i] < 0 " << q[i] << " coord: " << mesh.nodes[i].r << " " << mesh.nodes[i].z << endl;
		}*/
	}

	return 0;
}

double FEM_electro::V_in_point(double r, double z)
{
	double r1, r2, z1, z2, ksi_r, ksi_z;
	int node0 = 0, node_right = 2, node_top = 6;
	int cur_el = -1;
	for (int i = 0; i < mesh.kol_elements && cur_el == -1; i++)
	{
		r1 = mesh.nodes[mesh.elements[i].node_loc[node0]].r;
		r2 = mesh.nodes[mesh.elements[i].node_loc[node_right]].r;
		z1 = mesh.nodes[mesh.elements[i].node_loc[node0]].z;
		z2 = mesh.nodes[mesh.elements[i].node_loc[node_top]].z;
		
		if (r1 <= r && r <= r2 &&
			z1 <= z && z <= z2)
			cur_el = i;
	}
	if (cur_el == -1)
	{
		std::cout << "can't find (" << r << " " << z << ")\n";
		std::cout << "V is considered 0\n";
		return 0;
	}

	double V = 0;
	ksi_r = getKsi(r, r1, r2);
	ksi_z = getKsi(z, z1, z2);

	for (int i = 0; i < 9; i++)
	{
		V += basis[i%3](ksi_r) * basis[i/3](ksi_z) * q[mesh.elements[cur_el].node_loc[i]];
	}
	return V;
}

//int FEM::GetG_Loc(int el_id, double B) // получение локальной G
//{
//	double lambda = 1.0 / mesh.materials[mesh.elements[el_id].material].GetMu(IsLinear, B);
//	double hx = mesh.nodes[mesh.elements[el_id].node_loc[1]].x - mesh.nodes[mesh.elements[el_id].node_loc[0]].x,
//		   hy = mesh.nodes[mesh.elements[el_id].node_loc[3]].y - mesh.nodes[mesh.elements[el_id].node_loc[0]].y;
//
//	double koef1 = lambda * hy / 6.0 /hx ,
//		koef2 = lambda * hx / 6.0 / hy;
//	
//	A_loc[0][0] = koef1 * 2 + koef2 * 2;
//	A_loc[0][1] =-koef1 * 2 + koef2 * 1;
//	A_loc[0][2] = koef1 * 1 - koef2 * 2;
//	A_loc[0][3] =-koef1 * 1 - koef2 * 1;
//				  
//	A_loc[1][0] =-koef1 * 2 + koef2 * 1;
//	A_loc[1][1] = koef1 * 2 + koef2 * 2;
//	A_loc[1][2] =-koef1 * 1 - koef2 * 1;
//	A_loc[1][3] = koef1 * 1 - koef2 * 2;
//				 
//	A_loc[2][0] = koef1 * 1 - koef2 * 2;
//	A_loc[2][1] =-koef1 * 1 - koef2 * 1;
//	A_loc[2][2] = koef1 * 2 + koef2 * 2;
//	A_loc[2][3] =-koef1 * 2 + koef2 * 1;
//				 
//	A_loc[3][0] =-koef1 * 1 - koef2 * 1;
//	A_loc[3][1] = koef1 * 1 - koef2 * 2;
//	A_loc[3][2] =-koef1 * 2 + koef2 * 1;
//	A_loc[3][3] = koef1 * 2 + koef2 * 2;
//	return 0;
//}
//int FEM::GetG_Loc_test(int el_id) // получение локальной G
//{
//	double lambda = 1.0 / mu_test();
//	double 
//		hx = mesh.nodes[mesh.elements[el_id].node_loc[1]].x - mesh.nodes[mesh.elements[el_id].node_loc[0]].x,
//		hy = mesh.nodes[mesh.elements[el_id].node_loc[3]].y - mesh.nodes[mesh.elements[el_id].node_loc[0]].y;
//
//	double 
//		koef1 = lambda * hy / 6.0 / hx,
//		koef2 = lambda * hx / 6.0 / hy;
//
//	A_loc[0][0] = koef1 * 2 + koef2 * 2;
//	A_loc[0][1] = -koef1 * 2 + koef2 * 1;
//	A_loc[0][2] = koef1 * 1 - koef2 * 2;
//	A_loc[0][3] = -koef1 * 1 - koef2 * 1;
//
//	A_loc[1][0] = -koef1 * 2 + koef2 * 1;
//	A_loc[1][1] = koef1 * 2 + koef2 * 2;
//	A_loc[1][2] = -koef1 * 1 - koef2 * 1;
//	A_loc[1][3] = koef1 * 1 - koef2 * 2;
//
//	A_loc[2][0] = koef1 * 1 - koef2 * 2;
//	A_loc[2][1] = -koef1 * 1 - koef2 * 1;
//	A_loc[2][2] = koef1 * 2 + koef2 * 2;
//	A_loc[2][3] = -koef1 * 2 + koef2 * 1;
//
//	A_loc[3][0] = -koef1 * 1 - koef2 * 1;
//	A_loc[3][1] = koef1 * 1 - koef2 * 2;
//	A_loc[3][2] = -koef1 * 2 + koef2 * 1;
//	A_loc[3][3] = koef1 * 2 + koef2 * 2;
//	return 0;
//}
//int FEM::GetG_Loc_biquadratic(int el_id, double B) // получение локальной G
//{
//	double lambda = 1.0 / mesh.materials[mesh.elements[el_id].material].GetMu(IsLinear, B);
//	double
//		hx = mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[1]]].x - mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[0]]].x,
//		hy = mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[3]]].y - mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[0]]].y;
//
//	for (int i = 0; i < 9; i++)
//	{
//		for (int j = 0; j < 9; j++)
//		{
//			A_loc[i][j] = lambda * (hy / hx * G_quadratic[mu(i)][mu(j)] * M_quadratic[nu(i)][nu(j)]
//				+ hx / hy * M_quadratic[mu(i)][mu(j)] * G_quadratic[nu(i)][nu(j)]);
//		}
//	}
//	
//	return 0;
//}
//int FEM::GetG_Loc_biquadratic_test(int el_id) // получение локальной G
//{
//	double lambda = 1.0 / mu_test();
//	double
//		hx = mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[1]]].x - mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[0]]].x,
//		hy = mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[3]]].y - mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[0]]].y;
//
//	for (int i = 0; i < 9; i++)
//	{
//		for (int j = 0; j < 9; j++)
//		{
//			A_loc[i][j] = lambda * 
//				( hy / hx * G_quadratic[mu(i)][mu(j)] * M_quadratic[nu(i)][nu(j)]
//				+ hx / hy * M_quadratic[mu(i)][mu(j)] * G_quadratic[nu(i)][nu(j)]);
//		}
//	}
//
//	return 0;
//}
//
//int FEM::Getb_Loc(int el_id) // получение локального b
//{
//	double f = mesh.materials[mesh.elements[el_id].material].J;
//	double hx = mesh.nodes[mesh.elements[el_id].node_loc[1]].x - mesh.nodes[mesh.elements[el_id].node_loc[0]].x,
//		hy = mesh.nodes[mesh.elements[el_id].node_loc[3]].y - mesh.nodes[mesh.elements[el_id].node_loc[0]].y;
//	double coef = hx * hy * 9.0 / 36.0;
//	b_loc[0] = coef * f;
//	b_loc[1] = b_loc[0];
//	b_loc[2] = b_loc[0];
//	b_loc[3] = b_loc[0];
//	return 0;
//}
//int FEM::Getb_Loc_test(int el_id) // получение локального b
//{
//	double x1, x2, y1, y2;
//		x1 = mesh.nodes[mesh.elements[el_id].node_loc[0]].x;
//		x2 = mesh.nodes[mesh.elements[el_id].node_loc[1]].x;
//		y1 = mesh.nodes[mesh.elements[el_id].node_loc[0]].y;
//		y2 = mesh.nodes[mesh.elements[el_id].node_loc[3]].y;
//
//	double 
//		f1 = f_test(x1, y1),
//		f2 = f_test(x2, y1),
//		f3 = f_test(x1, y2),
//		f4 = f_test(x2,y2);
//	double 
//		hx = x2-x1,
//		hy = y2-y1;
//	double coef = hx * hy  / 36.0;
//	b_loc[0] = 4*f1 + 2*f2 + 2*f3 + f4;
//	b_loc[1] = 2*f1 + 4*f2 + f3 + 2*f4;
//	b_loc[2] = 2*f1 + f2 + 4*f3 + 2*f4;
//	b_loc[3] = f1 + 2*f2 + 2*f3 + 4*f4;
//	return 0;
//}
//int FEM::Getb_Loc_biquadratic(int el_id) // получение локального b
//{
//	double f = mesh.materials[mesh.elements[el_id].material].J;
//	double 
//		hx = mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[1]]].x - mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[0]]].x,
//		hy = mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[3]]].y - mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[0]]].y;
//	double coef = f * hx * hy;
//	for (int i = 0; i < 9; i++)
//	{
//		b_loc[i] = 0;
//		for (int j = 0; j < 9; j++)
//		{
//			b_loc[i] += coef * M_quadratic[mu(i)][mu(j)] * M_quadratic[nu(i)][nu(j)];
//		}
//	}
//
//	return 0;
//}
//int FEM::Getb_Loc_biquadratic_test(int el_id) // получение локального b
//{
//	vector<double> f(9);
//	for (int i = 0; i < 9; i++)
//	{
//		f[i] = f_test(mesh.nodes[mesh.elements[el_id].node_loc[i]].x, mesh.nodes[mesh.elements[el_id].node_loc[i]].y);
//	}
//	double
//		hx = mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[1]]].x - mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[0]]].x,
//		hy = mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[3]]].y - mesh.nodes[mesh.elements[el_id].node_loc[conner_loc[0]]].y;
//	double coef = hx * hy;
//	for (int i = 0; i < 9; i++)
//	{
//		b_loc[i] = 0;
//		for (int j = 0; j < 9; j++)
//		{
//			b_loc[i] += f[j] * coef * M_quadratic[mu(i)][mu(j)] * M_quadratic[nu(i)][nu(j)];
//		}
//	}
//
//	return 0;
//}
//
//int FEM::Init(string path)
//{
//	conner_loc = {0, 1, 2, 3};
//	basis = { linear1, linear2 };
//	//mesh.MakeLinear_test(path);
//	IsLinear = 1;
//	IsLinear = mesh.ReadMesh(path);
//	slau.GeneratePortret(mesh);
//	q.resize(mesh.kol_nodes);
//	b_loc.resize(4);
//	A_loc.resize(4);
//	for (int i = 0; i < 4; i++)
//	{
//		A_loc[i].resize(4);
//	}
//	return 0;
//}
//int FEM::Init_biquadratic(string path)
//{
//	conner_loc = { 0, 2, 6, 8 };
//	basis = { quadratic1, quadratic2, quadratic3};
//	IsLinear = mesh.ReadMesh(path);
//	mesh.MakeBiquadratic(path);
//	slau.GeneratePortret(mesh);
//	q.resize(mesh.kol_nodes);
//	b_loc.resize(9);
//	A_loc.resize(9);
//	for (int i = 0; i < 9; i++)
//	{
//		A_loc[i].resize(9);
//	}
//	return 0;
//}
//
//int FEM::SolveTask(string path)
//{
//	Init(path);
//
//	for (int i = 0; i < mesh.kol_elements; i++)
//	{
//		Getb_Loc(i);
//		GetG_Loc(i);
//		slau.AddLocal(A_loc, b_loc, i, mesh);
//	}
//	slau.SetS1(mesh);
//
//	slau.ShowMatrix();
//
//	cgm_solver solver1(slau.ia, slau.ja, slau.di, slau.a, slau.b);
//	solver1.llt_preconditioning(q);
//	solver1.clear_all();
//
//	if (!IsLinear)
//	{
//		q_prev.resize(q.size());
//		Node centre;
//		double B, _, discrepancy;
//		pair<double, double> B_vec;
//		int stop = 0;
//		for (int i = 0; i < MAX_ITER && !stop; i++)
//		{
//			slau.ClearValues();
//			for (int i = 0; i < mesh.kol_elements; i++)
//			{
//				centre = mesh.GetCentre(i);
//				B_vec = GetBInPoint(centre.x, centre.y, _, i);
//				B = sqrt(B_vec.first * B_vec.first + B_vec.second * B_vec.second);
//
//				Getb_Loc(i);
//				GetG_Loc(i, B);
//				slau.AddLocal(A_loc, b_loc, i, mesh);
//			}
//			slau.SetS1(mesh);
//
//			if (i != 0)
//			{
//				discrepancy = sqrt(CheckRes());
//				cout << "--------------- discrepancy = " << discrepancy << " ---------------" << endl << endl;
//				if (discrepancy < EPS_NONLINEAR)
//					stop = 1;
//			}
//			if (!stop)
//			{
//				cout << "--------------- NONLINEAR ITERATION " << i << " ---------------" << endl;
//				cgm_solver solver(slau.ia, slau.ja, slau.di, slau.a, slau.b);
//				q_prev.swap(q);
//				solver.llt_preconditioning(q);
//				for (int i = 0; i < q.size(); i++)
//				{
//					q[i] = w * q[i] + (1 - w) * q_prev[i];
//				}
//			}
//		}
//	}
//
//	return 0;
//}
//int FEM::SolveTask_test(string path)
//{
//	Init(path);
//
//	for (int i = 0; i < mesh.kol_elements; i++)
//	{
//		Getb_Loc_test(i);
//		GetG_Loc_test(i);
//		slau.AddLocal(A_loc, b_loc, i, mesh);
//	}
//	slau.SetS1_test(mesh);
//
//	slau.ShowMatrix();
//
//	cgm_solver solver(slau.ia, slau.ja, slau.di, slau.a, slau.b);
//	solver.llt_preconditioning(q);
//	return 0;
//}
//int FEM::SolveTask_biquadratic(string path)
//{
//	Init_biquadratic(path);
//	cout << "Build mesh\n";
//	for (int i = 0; i < mesh.kol_elements; i++)
//	{
//		Getb_Loc_biquadratic(i);
//		GetG_Loc_biquadratic(i);
//		slau.AddLocal(A_loc, b_loc, i, mesh);
//	}
//	slau.SetS1(mesh);
//	cout << "Start solve SLAU\n";
//	cgm_solver solver1(slau.ia, slau.ja, slau.di, slau.a, slau.b);
//	solver1.llt_preconditioning(q);
//	solver1.clear_all();
//
//	if (!IsLinear)
//	{
//		q_prev.resize(q.size());
//		Node centre;
//		double B, _, discrepancy = 1e30, prev_discrepancy = 1e30;
//		pair<double, double> B_vec;
//		int stop = 0;
//		for (int i = 0; i < MAX_ITER && !stop; i++)
//		{
//			slau.ClearValues();
//			for (int i = 0; i < mesh.kol_elements; i++)
//			{
//				centre = mesh.GetCentre(i);
//				B_vec = GetBInPoint(centre.x, centre.y, _, i);
//				B = sqrt(B_vec.first * B_vec.first + B_vec.second * B_vec.second);
//
//				Getb_Loc_biquadratic(i);
//				GetG_Loc_biquadratic(i, B);
//				slau.AddLocal(A_loc, b_loc, i, mesh);
//			}
//			slau.SetS1(mesh);
//
//			if (i != 0)
//			{
//				prev_discrepancy = discrepancy;
//				discrepancy = sqrt(CheckRes());
//				cout << "--------------- discrepancy = " << discrepancy << " ---------------" << endl << endl;
//				if (discrepancy < EPS_NONLINEAR)
//					stop = 1;
//				if (prev_discrepancy < discrepancy && w > 1e-5)
//					w /= 2;
//			}
//			if (!stop)
//			{
//				cout << "--------------- NONLINEAR ITERATION " << i << " ---------------" << endl;
//				cout << "Start solve SLAU\n";
//				cgm_solver solver(slau.ia, slau.ja, slau.di, slau.a, slau.b, 1e-15);
//				q_prev.swap(q);
//				solver.llt_preconditioning(q);
//				//solver.diag_preconditioning(q);
//
//				for (int i = 0; i < q.size(); i++)
//				{
//					q[i] = w * q[i] + (1 - w) * q_prev[i];
//				}
//			}
//		}
//	}
//
//	return 0;
//}
//int FEM::SolveTask_biquadratic_test(string path)
//{
//	Init_biquadratic(path);
//
//	for (int i = 0; i < mesh.kol_elements; i++)
//	{
//		Getb_Loc_biquadratic_test(i);
//		GetG_Loc_biquadratic_test(i);
//		slau.AddLocal(A_loc, b_loc, i, mesh);
//		slau.AddLocal(A_loc, b_loc, i, mesh);
//	}
//	slau.SetS1_test(mesh);
//
//	//Solver solver(slau.ia, slau.ja, slau.di, slau.a, slau.a, slau.b);
//	//solver.BSG();
//	//solver.getx0(q);
//
//	cgm_solver solver(slau.ia, slau.ja, slau.di, slau.a, slau.b);
//	solver.llt_preconditioning(q);
//	return 0;
//}
