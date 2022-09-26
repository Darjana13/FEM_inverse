#include "FEM.h"

int SLAU::GeneratePortret(Mesh& mesh)
{
	int N = mesh.kol_nodes, Kel = mesh.kol_elements;
	int N_loc = mesh.elements[0].node_loc.size();
	ia.resize(N + 1);
	
	vector<set<int>> map(N);
	int k = 0;
	for (auto elem : mesh.elements)
	{
		for (auto i : elem.node_loc)
			for (auto j : elem.node_loc)
				if (i > j)
					map[i].insert(j);
	}
	ia[0] = 0;
	for (int i = 0; i < N; i++)
	{
		ia[i + 1] = ia[i] + map[i].size();
	}

	ja.resize(ia[N]);
	for (int i = 0; i < N; i++)
	{
		set <int> ::iterator it = map[i].begin();
		for (int j = 0; it != map[i].end(); it++, j++)
		{
			ja[ia[i] + j] = *(it);
		}
	}
	
	au.resize(ja.size());
	al.resize(ja.size());

	di.resize(N);
	b.resize(N);

	return 0;
}

void SLAU::ClearValues()
{
	fill(au.begin(), au.end(), 0.0);
	fill(al.begin(), al.end(), 0.0);

	fill(b.begin(), b.end(), 0.0);
	fill(di.begin(), di.end(), 0.0);
}

int SLAU::ShowMatrix()
{
	ofstream ofp("Matrix before solution.txt");
	int N = di.size();
	ofp << "A: row, col - value" << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			int temp = ia[i];
			for (int k = temp; k < ia[i + 1]; k++)
			{
				if (ja[k] == j)
				{
					ofp << "   " << i << ", " << ja[k] << " - au " << au[k] << " al " << al[k] << "\n";
					break;
				}
				
			}
		}
	}
	for (int i = 0; i < N; i++)
	{
		
		ofp << "b: " << i << " - " << b[i] << "\n";
		
	}
	for (int i = 0; i < N; i++)
	{

		ofp << "di: " << i << " - " << di[i] << "\n";

	}
	return 0;
}

int SLAU::AddLocal(vector<vector<double>>& A_loc, vector<double>& b_loc, int el_id, Mesh& mesh)
// внесение локальных A, b  в глобальную СЛАУ
{
	vector<int> L = mesh.elements[el_id].node_loc;
	int n_loc = mesh.elements[el_id].node_loc.size(); // размерность локальной матрицы
	
	for (int i = 0; i < n_loc; i++)
	{
		di[L[i]] += A_loc[i][i];
		b[L[i]] += b_loc[i];
	}
	
	for (int i = 0; i < n_loc; i++)
	{
		int temp = ia[L[i]];
		for (int j = 0; j < i; j++)
		{
			for (int k = temp; k < ia[L[i] + 1]; k++)
			{
				if (ja[k] == L[j])
				{
					au[k] += A_loc[j][i];	
					al[k] += A_loc[i][j];
					break;
				}
			}
		}
	}

	return 0;
}

double GetUTest(double r, double z)
{
	return z*z*z + r * r * r;
}

int SLAU::SetS1_null(Mesh& mesh) // учет первых краевых
{
	for (int i = 0; i < mesh.kol_S1nodes; i++)
	{
		di[mesh.S1[i]] = 1;
		b[mesh.S1[i]] = 0;
		
		//cout << " i " << i << " r " << mesh.nodes[mesh.S1[i]].r << mesh.nodes[mesh.S1[i]].z << endl;

		int node_id = mesh.S1[i];
		for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
		{
			al[k] = 0;
			au[k] = 0;
		}
		for (int k = 0; k < ja.size(); k++)
		{
			if (ja[k] == node_id)
			{
				al[k] = 0;
				au[k] = 0;
			}
		}
	}
	return 0;
}


int SLAU::SetS1(Mesh& mesh) // учет первых краевых
{
	for (int i = 0; i < mesh.kol_S1nodes; i++)
	{
		di[mesh.S1[i]] = 1;
		b[mesh.S1[i]] = 0;
		//b[mesh.S1[i]] = GetUTest(mesh.nodes[mesh.S1[i]].r, mesh.nodes[mesh.S1[i]].z);
		
		
		int node_id = mesh.S1[i];
		for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
		{
			al[k] = 0;
		}
		for (int k = 0; k < ja.size(); k++)
		{
			if (ja[k] == node_id)
			{
				au[k] = 0;
			}
		}
	}
	return 0;
}
//int SLAU::SetS1_test(Mesh& mesh) // учет первых краевых
//{
//	int NS1 = mesh.kol_S1nodes;
//	double max = 0;
//	for (int i = 0; i < di.size(); i++)
//	{
//		if (di[i] > max)
//			max = di[i];
//	}
//	for (int i = 0; i < NS1; i++)
//	{
//		di[mesh.S1[i]] = max * BIG_VALUE;
//		b[mesh.S1[i]] = max * BIG_VALUE * u_test(mesh.nodes[mesh.S1[i]].x, mesh.nodes[mesh.S1[i]].y);
//		cout << mesh.nodes[mesh.S1[i]].x << " " << mesh.nodes[mesh.S1[i]].y << " " << u_test(mesh.nodes[mesh.S1[i]].x, mesh.nodes[mesh.S1[i]].y) << endl;
//	}
//	return 0;
//}
