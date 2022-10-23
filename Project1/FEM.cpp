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
		in >> receivers[i].r[0] >> receivers[i].z[0] >> receivers[i].r[1] >> receivers[i].z[1] >> receivers[i].V_true;
	}
	in.close();

	in.open(path + "VEL.txt");
	in >> kol_vel;
	VELs.resize(kol_vel);
	int cur_el[2];
	double r1, r2, z1, z2;
	for (int i = 0; i < kol_vel; i++)
	{
		in >> VELs[i].r[0] >> VELs[i].z[0] >> VELs[i].r[1] >> VELs[i].z[1] >> VELs[i].J[0];
		VELs[i].J[0] /= (2.0 * M_PI);
		VELs[i].J[1] = -VELs[i].J[0];
		for (int k = 0; k < 2; k++)
		//for (int k = 0; k < 1; k++)
		{
			cur_el[k] = FindElem(VELs[i].r[k], VELs[i].z[k]);
			
			double hr = mesh.nodes[mesh.elements[cur_el[k]].node_loc[2]].r - mesh.nodes[mesh.elements[cur_el[k]].node_loc[0]].r;
			double hz = mesh.nodes[mesh.elements[cur_el[k]].node_loc[6]].z - mesh.nodes[mesh.elements[cur_el[k]].node_loc[0]].z;

			for (int j = 0; j < 9; j++)
			{
				f_ist[mesh.elements[cur_el[k]].node_loc[j]] = VELs[i].J[k] / hr/hr/hz;
			}
			/*f_ist[mesh.elements[cur_el].node_loc[3]] = VELs[i].J[k] / hr / hr / hz;
			f_ist[mesh.elements[cur_el].node_loc[4]] = VELs[i].J[k] / hr / hr / hz;
			f_ist[mesh.elements[cur_el].node_loc[6]] = VELs[i].J[k] / hr / hr / hz;
			f_ist[mesh.elements[cur_el].node_loc[7]] = VELs[i].J[k] / hr / hr / hz;*/
		}
		if (cur_el[0] == cur_el[1])
		{
			cout << "Error: VEL is placed in one elem " << endl;
			exit(1);
		}
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
			{
				cout << "nan" << endl;
				exit(1);
			}
		}
	}
	
	return 0;
}

double FEM_electro::GetUTest(double r, double z)
{
	//return z*z*z+r*r*r;
	return r == 0 ? -VELs[0].J[0] / mesh.sigma[0] / (1e-2) : -VELs[0].J[0] / (r * mesh.sigma[0]);
}

double FEM_electro::GetF(int node_id)
{
	if (f_ist.count(node_id))
		return f_ist[node_id];
	return 0;
}

double FEM_electro::GetF(double r, double z) // возвращает только в узле J
{
	//return r == 0 ? -VELs[0].J[0] / (1e-6) : -VELs[0].J[0] / (r * r * r);
	//return -6*z-9*r;

	for (int i = 0; i < kol_vel; i++)
	{
		if (VELs[i].r[0] == r && VELs[i].z[0] == z)
		{
			cout << r << " " << z << " " << VELs[i].J[0] << endl;
			return VELs[i].J[0];
		}
		if (VELs[i].r[1] == r && VELs[i].z[1] == z)
			return VELs[i].J[1];
	}
	return 0;
}

int FEM_electro::Getb_Loc_biquadratic(int el_id) // получение локального b
{
	vector<double> f(9);
	bool not_null = 0;
	double hr = 0, hz = 0;
	double r_cur;

	for (int i = 0; i < 9; i++)
	{
		//f[i] = GetF(mesh.nodes[mesh.elements[el_id].node_loc[i]].r, mesh.nodes[mesh.elements[el_id].node_loc[i]].z);
		f[i] = GetF(mesh.elements[el_id].node_loc[i]);
		if (f[i] != 0)
		{
			// объемный источник, J = I/2*pi
			hr = mesh.nodes[mesh.elements[el_id].node_loc[2]].r - mesh.nodes[mesh.elements[el_id].node_loc[0]].r;
			hz = mesh.nodes[mesh.elements[el_id].node_loc[6]].z - mesh.nodes[mesh.elements[el_id].node_loc[0]].z;
			not_null = 1;
		}
	}

	if (not_null == 1)
	{
		r_cur = mesh.nodes[mesh.elements[el_id].node_loc[0]].r;
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

int FEM_electro::DirectTask(int mat_id, double h_sigma)
{
	double save_sigma = mesh.sigma[mat_id];
	mesh.sigma[mat_id] = h_sigma;
	DirectTask();
	mesh.sigma[mat_id] = save_sigma;
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

		if(i % 10000 == 0)
			cout << "el " << i << " from " << mesh.kol_elements << endl;
	}
	cout << "end make matrix " << endl;

	cout << "start Set S1 " << endl;
	//slau.SetS1(mesh);
	slau.SetS1_null(mesh);
	cout << "end Set S1 " << endl;

	cgm_solver solver(slau.ia, slau.ja, slau.di, slau.au, slau.b, 1e-10);
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

	/*vector<double> q_ist(mesh.kol_nodes);
	for (int i = 0; i < mesh.kol_nodes; i++)
	{
		q_ist[i] = VELs[0].J[0] / mesh.sigma[0] / (sqrt(mesh.nodes[i].r * mesh.nodes[i].r + mesh.nodes[i].z * mesh.nodes[i].z));
	}
	solver.A.Ax(q_ist, b_ist);
	cout << " A*q_ist ";
	for (int i = 0; i < mesh.kol_nodes; i++)
	{
		if (abs(slau.b[i] - b_ist[i]) > 1e-5)
			cout << " i " << i << " q_ist " << q_ist[i] << " q " << q[i] << " b[i] " << slau.b[i] << " Aq[i] " << b_ist[i] << endl;
	}*/

	/*double V_start, V_end;
	for (int i = 0; i < kol_rec; i++)
	{
		V_start = V_in_point(receivers[i].r[0], receivers[i].z[0]);
		V_end = V_in_point(receivers[i].r[1], receivers[i].z[1]);
		receivers[i].V = V_start - V_end;
	}*/

	return 0;
}

int FEM_electro::V_in_rec(vector<double> &V)
{
	V.resize(kol_rec);
	double V_start, V_end;
	for(int i = 0; i < kol_rec; i++)
	{
		 V_start = V_in_point(receivers[i].r[0], receivers[i].z[0]);
		 V_end = V_in_point(receivers[i].r[1], receivers[i].z[1]);

		 V[i] = V_start - V_end;
	}
	return 0;
}

int FEM_electro::PrintLine(string filename, double z, double r)
{
	ofstream line_out(filename);
	for (int i = 0; i < mesh.kol_nodes; i += 1)
	{
		if (mesh.nodes[i].r > r && mesh.nodes[i].z == z)
		{
			line_out << mesh.nodes[i].r << "\t" << q[i] << endl;
		}
	}
	return 0;
}

int FEM_electro::PrintField(string filename)
{
	ofstream out(filename);
	for (int i = 0; i < mesh.kol_nodes; i+=1)
	{
		out << mesh.nodes[i].r << " " << mesh.nodes[i].z << " " << q[i] << endl;
	}
	return 0;
}

int FEM_electro::FindElem(double r, double z)
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
	return cur_el;
}

double FEM_electro::V_in_point(double r, double z)
{
	double r1, r2, z1, z2, ksi_r, ksi_z;
	int cur_el = FindElem(r, z);
	if (cur_el == -1)
	{
		std::cout << "can't find (" << r << " " << z << ")\n";
		std::cout << "V is considered 0\n";
		return 0;
	}
	r1 = mesh.nodes[mesh.elements[cur_el].node_loc[0]].r;
	r2 = mesh.nodes[mesh.elements[cur_el].node_loc[2]].r;
	z1 = mesh.nodes[mesh.elements[cur_el].node_loc[0]].z;
	z2 = mesh.nodes[mesh.elements[cur_el].node_loc[6]].z;

	double V = 0;
	ksi_r = getKsi(r, r1, r2);
	ksi_z = getKsi(z, z1, z2);

	for (int i = 0; i < 9; i++)
	{
		V += basis[i%3](ksi_r) * basis[i/3](ksi_z) * q[mesh.elements[cur_el].node_loc[i]];
	}
	return V;
}
