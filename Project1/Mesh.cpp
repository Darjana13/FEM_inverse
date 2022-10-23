#include "FEM.h"

int PointsOnAxis(std::ifstream& in, std::vector<double>& all_X, int count_x) // qudratic
{
	double X, kx;
	int Nx;
	for (int curr_count_x = 0; curr_count_x < count_x - 1; )
	{
		in >> X >> Nx >> kx;
		double hx;
		if (kx == 1)
		{
			hx = (X - all_X[curr_count_x]) / Nx;
			for (int p = 1; p < Nx; p++)
			{
				all_X[curr_count_x + 2 * p] = all_X[curr_count_x] + hx * p;
				all_X[curr_count_x + 2 * p - 1] = (all_X[curr_count_x + 2 * p - 2] + all_X[curr_count_x + 2 * p]) / 2;
			}
			curr_count_x += 2*Nx;
		}
		else
		{
			hx = (X - all_X[curr_count_x]) * (kx - 1) / (pow(kx, Nx) - 1);
			curr_count_x++;
			for (int p = 0; p < Nx - 1; curr_count_x += 2, p++)
			{
				all_X[curr_count_x + 1] = all_X[curr_count_x - 1] + hx * pow(kx, p);
				all_X[curr_count_x] = (all_X[curr_count_x + 1] + all_X[curr_count_x - 1]) / 2;
			}
			curr_count_x+=1;

		}
		all_X[curr_count_x] = X;
		all_X[curr_count_x - 1] = (all_X[curr_count_x] + all_X[curr_count_x - 2]) / 2;
	}
	return 0;
}

int Mesh::ReadMat(string path)
{
	ifstream in(path + "mat.txt");
	ofstream out;
	in >> kol_mat;
	sigma.resize(kol_mat);
	vector<double> h_layer(kol_mat);
	for (int i = 0; i < kol_mat; i++)
	{
		in >> h_layer[i] >> sigma[i];
	}
	in.close();
	int mat_id = 0;
	for (int i = 0; i < kol_elements; i++)
	{
		if (GetCentre(i).z > h_layer[mat_id])
			mat_id++;
		elements[i].material = mat_id;
	}

	out.open(path + "elem.txt");
	out << kol_elements << endl;
	for (int i = 0; i < kol_elements; i++)
	{
		for (int j = 0; j < elements[i].node_loc.size(); j++)
		{
			out << elements[i].node_loc[j] << " ";
		}
		out << elements[i].material << endl;
	}
	out.close();
	return 0;
}

int Mesh::ReadMesh(string path)
{
	CreateMesh(path);
	ReadMat(path);
	return 0;
}

int Mesh::CreateMesh(string path)
{
	std::ofstream out;
	out.precision(15);
	std::vector<double> all_R, all_Z;
	std::ifstream in(path + "grid.txt");
	double R, Z, kr, kz;
	int Nr, Nz;
	int count_r, count_z;
	vector<double> h_layer, sigma;
	in >> count_r >> count_z;
	all_R.resize(count_r*2-1);
	all_Z.resize(count_z*2-1);

	in >> all_R[0] >> all_Z[0];
	PointsOnAxis(in, all_R, count_r * 2 - 1);
	PointsOnAxis(in, all_Z, count_z * 2 - 1);
	in.close();

	kol_nodes = (count_z * 2 - 1) * (count_r * 2 - 1);
	nodes.resize(kol_nodes);
	for (int i = 0; i < count_z * 2 - 1; i++)
	{
		for (int k = 0; k < count_r * 2 - 1; k++)
		{
			nodes[i * (count_r * 2 - 1) + k].r = all_R[k];
			nodes[i * (count_r * 2 - 1) + k].z = all_Z[i];
		}
	}

	out.open(path + "xyz.txt");
	out << kol_nodes << endl;
	for (int i = 0; i < count_z * 2 - 1; i++)
	{
		for (int k = 0; k < count_r * 2 - 1; k++)
		{
			out << all_R[k] << "\t" << all_Z[i] << "\n";
		}
	}
	out.close();

	out.open(path + "r.txt");
		for (int k = 0; k < count_r * 2 - 1; k++)
		{
			out << all_R[k] << "\n";
		}
	out.close();

	out.open(path + "z.txt");
	for (int i = 0; i < count_z * 2 - 1; i++)
	{
			out << all_Z[i] << "\n";
	}
	out.close();

	int el_id = 0;
	kol_elements = (count_z - 1) * (count_r - 1);
	elements.resize(kol_elements);
	for (int i = 0; i < 2 * count_z - 2; i += 2)
	{
		for (int j = 0; j < 2 * count_r - 2; j += 2)
		{
			elements[el_id].node_loc =
			{ i * (2 * count_r - 1) + j,      i * (2 * count_r - 1) + j + 1,      i * (2 * count_r - 1) + j + 2,
			  (i + 1) * (2 * count_r - 1) + j, (i + 1) * (2 * count_r - 1) + j + 1, (i + 1) * (2 * count_r - 1) + j + 2,
			  (i + 2) * (2 * count_r - 1) + j, (i + 2) * (2 * count_r - 1) + j + 1, (i + 2) * (2 * count_r - 1) + j + 2,
			};
			el_id++;
		}
	}
	

	kol_S1nodes = (2 * count_z - 1 + 2 * count_r - 1);
	S1.resize(kol_S1nodes);
	int S1_id = 0;
	// left
	//for (int i = 0; i < 2 * count_z - 1; i++)
	//{
	//	S1[S1_id] = i * (2 * count_r - 1);
	//	S1_id++;
	//}
	// top
	//for (int i = 1; i < 2 * count_r; i++)
	//{
	//	S1[S1_id] = (2 * count_z - 1) * (2 * count_r - 1) - i;
	//	S1_id++;
	//}
	// right
	for (int i = 1; i < 2 * count_z - 1 + 1; i++)
	{
		S1[S1_id] = i * (2 * count_r - 1) - 1;
		S1_id++;
	}
	// bottom
	for (int i = 0; i < 2 * count_r - 1; i++)
	{
		S1[S1_id] = i;
		S1_id++;
	}

	out.open(path + "S1.txt");
	out << kol_S1nodes << endl;
	for (int i = 0; i < kol_S1nodes; i++)
	{
		out << S1[i] << endl;
	}
	out.close();

	return 0;
}

Node Mesh::GetCentre(int el_id)
{
	double x = 0, y = 0;
	int i;
	for (i = 0; i < elements[el_id].node_loc.size(); i++)
	{
		x += nodes[elements[el_id].node_loc[i]].r;
		y += nodes[elements[el_id].node_loc[i]].z;
	}
	x /= i;
	y /= i;
	return Node(x, y);
}