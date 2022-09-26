#include "FEM.h"

int PointsOnAxis(std::ifstream& in, std::vector<double>& all_X, double count_x) // qudratic
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
			curr_count_x+=3;

			/*hx = 0.01;
			kx = 1.1;
			curr_count_x++;
			for (int p = 0; p < Nx - 1; curr_count_x += 2, p++)
			{
				all_X[curr_count_x + 1] = all_X[curr_count_x - 1] + hx * pow(kx, p);
				all_X[curr_count_x] = (all_X[curr_count_x + 1] + all_X[curr_count_x - 1]) / 2;
				if (all_X[curr_count_x + 1] >= 100)
					int iiiii = 0;
			}
			curr_count_x += 3;*/
		}
		all_X[curr_count_x] = X;
		all_X[curr_count_x - 1] = (all_X[curr_count_x] + all_X[curr_count_x - 2]) / 2;
	}
	return 0;
}


//int Mesh::ReadMesh(string path)
//{
//	string tmpstr;
//	int tmpint;
//	long tmplong;
//	double tmpdbl;
//
//	int IsLinear = 1;
//
//	/*inf2tr.dat(текстовый файл) :
//		kuzlov - число узлов сетки
//		ktr - число конечных элементов
//		kt1 - число узлов с первыми краевыми условиями
//		(остальные параметры для выполнения лабораторной работы не используются)
//		islau = 0 indku1 = 0 indfro = 1
//		kuzlov = 1080    ktr = 1014    kt1 = 94  kreb2 = 0  kreb3 = 0
//		kisrr1 = 2 kisrr2 = 2 kisrr3 = 2  kbrsr = 8
//		kreb4 = 0*/
//
//	string file = path + "inf2tr.dat";
//	ifstream inf(file);
//	if (!inf)
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	inf >> tmpstr >> tmpstr >> tmpstr >> tmpstr >> tmpstr >> tmpstr >> tmpstr;
//	inf >> kol_nodes >> tmpstr >> kol_elements >> tmpstr >> kol_S1nodes;
//	inf.close();
//	inf.clear();
//	elements.resize(kol_elements);
//	nodes.resize(kol_nodes);
//	S1.resize(kol_S1nodes);
//
//	/*nvtr.dat(двоичный файл) :
//		число записей - ktr
//		структура i - й записи : 6 * long(i1, i2, i3, i4, 0, 1),
//		где i1, i2, i3, i4 - номера(глобальные) четырех вершин i - го прямоугольника*/
//	file = path + "nvtr.dat";
//	inf.open(file, ios::binary);
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	for (int i = 0; i < kol_elements; i++)
//	{
//		// потому что в файле записаны сначала верхние точки элемента, потом нижние
//		inf.read((char*)&elements[i].node_loc[2], sizeof(long)); if (inf.gcount() != sizeof(long)) break;
//		inf.read((char*)&elements[i].node_loc[3], sizeof(long)); if (inf.gcount() != sizeof(long)) break;
//		inf.read((char*)&elements[i].node_loc[0], sizeof(long)); if (inf.gcount() != sizeof(long)) break;
//		inf.read((char*)&elements[i].node_loc[1], sizeof(long)); if (inf.gcount() != sizeof(long)) break;
//		inf.read((char*)&tmplong, sizeof(long)); if (inf.gcount() != sizeof(long)) break;
//		inf.read((char*)&tmplong, sizeof(long)); if (inf.gcount() != sizeof(long)) break;
//
//		elements[i].node_loc[0]--;
//		elements[i].node_loc[1]--;
//		elements[i].node_loc[2]--;
//		elements[i].node_loc[3]--;
//	}
//	inf.close();
//	inf.clear();
//
//	/*nvkat2d.dat(двоичный файл) :
//		число записей - ktr
//		структура i - й записи : 1 * long(номер материала i - го прямоугольника в соответствии с табл.1)*/
//	file = path + "nvkat2d.dat";
//	inf.open(file, ios::binary);
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	kol_mat = 0;
//	for (int i = 0; i < kol_elements; i++)
//	{
//		inf.read((char*)&elements[i].material, sizeof(long)); if (inf.gcount() != sizeof(long)) break;
//		if (elements[i].material > kol_mat)
//			kol_mat = elements[i].material;
//		elements[i].material--;
//	}
//	inf.close();
//	inf.clear();
//
//	/*rz.dat(двоичный файл) :
//		число записей - kuzlov
//		структура i - й записи : 2 * double(x, y)
//		где x, y - (x, y) - координаты i - й вершины*/
//	file = path + "rz.dat";
//	inf.open(file, ios::binary);
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	for (int i = 0; i < kol_nodes; i++)
//	{
//		inf.read((char*)&nodes[i].x, sizeof(double)); if (inf.gcount() != sizeof(double)) break;
//		inf.read((char*)&nodes[i].y, sizeof(double)); if (inf.gcount() != sizeof(double)) break;
//	}
//	inf.close();
//	inf.clear();
//
//	/*l1.dat(двоичный файл) :
//		число записей - kt1
//		структура i - й записи : 1 * long(номер i - й вершины с первым нулевыми краевым условием)*/
//	file = path + "l1.dat";
//	inf.open(file, ios::binary);
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	for (int i = 0; i < kol_S1nodes; i++)
//	{
//		inf.read((char*)&S1[i], sizeof(long)); if (inf.gcount() != sizeof(long)) break;
//		S1[i]--;
//	}
//	inf.close();
//	inf.clear();
//
//	/*mu_test(текстовый файл)
//		число строк - число материалов
//		структура i - й строки : номер_материала  соответствующее_значение_мю*/
//	materials.resize(kol_mat);
//	file = path + "mu";
//	inf.open(file);
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	ifstream exist;
//	for (int i = 0; i < kol_mat; i++)
//	{
//		inf >> tmpint >> tmpdbl;
//		file = path + "mu.00" + NumberToString(tmpint);
//		exist.open(file);
//		if (!exist)
//			materials[tmpint - 1].mu = MU0 * tmpdbl;
//		else
//		{
//			materials[tmpint - 1].mu = MU0 * tmpdbl;
//			exist.close();
//			exist.clear();
//			IsLinear = 0;
//			materials[tmpint - 1].SetMu(file);
//		}
//	}
//	inf.close();
//	inf.clear();
//
//	/*toku(текстовый файл)
//		число строк - число материалов
//		структура i - й строки : номер_материала  соответствующее_значение_тока*/
//	file = path + "toku";
//	inf.open(file);
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	for (int i = 0; i < kol_mat; i++)
//	{
//		inf >> tmpint >> tmpdbl;
//		materials[tmpint - 1].J = tmpdbl;
//	}
//	inf.close();
//	inf.clear();
//
//	return IsLinear;
//}

//int Mesh::MakeBiquadratic(string path)
//{
//	string file = path + "r.dat";
//	ifstream inf;
//	double tmpdbl;
//	int tmpint = 0, x_kol, y_kol;
//	vector<double> x, y;
//
//	inf.open("rasm");
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	inf >> x_kol >> y_kol;
//	inf.close();
//	inf.clear();
//
//	x.resize(2*x_kol - 1);
//	inf.open(file, ios::binary);
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	inf.read((char*)&x[0], sizeof(double));
//	for (int i = 1; i < 2*x_kol - 2; i+=2 )
//	{
//		inf.read((char*)&tmpdbl, sizeof(double)); if (inf.gcount() != sizeof(double)) break;
//		x[i] = (tmpdbl + x[i - 1]) / 2.0;
//		x[i + 1] = tmpdbl;
//	}
//	inf.close();
//	inf.clear();
//
//	y.resize(2 * y_kol - 1);
//	file = path + "z.dat";
//	inf.open(file, ios::binary);
//	if (!inf.is_open())
//	{
//		cout << "error Can't open " << file << endl;
//		exit(1);
//	}
//	inf.read((char*)&y[0], sizeof(double));
//	for (int i = 1; i < 2 * y_kol - 2; i += 2)
//	{
//		inf.read((char*)&tmpdbl, sizeof(double)); if (inf.gcount() != sizeof(double)) break;
//		y[i] = (tmpdbl + y[i - 1]) / 2.0;
//		y[i + 1] = tmpdbl;
//	}
//	inf.close();
//	inf.clear();
//
//	int node_id = 0;
//	kol_nodes = (2 * y_kol - 1) * (2 * x_kol - 1);
//	nodes.resize(kol_nodes);
//
//	for (int i = 0; i < 2*y_kol - 1; i++)
//	{
//		for (int j = 0; j < 2*x_kol - 1; j++)
//		{
//			nodes[node_id].x = x[j];
//			nodes[node_id].y = y[i];
//			node_id++;
//		}
//	}
//
//	int el_id = 0;
//	for (int i = 0; i < 2 * y_kol - 2; i+=2)
//	{
//		for (int j = 0; j < 2 * x_kol - 2; j+=2)
//		{
//			
//			elements[el_id].node_loc = 
//			{ i       * (2 * x_kol - 1) + j,      i * (2 * x_kol - 1) + j + 1,      i * (2 * x_kol - 1) + j + 2,
//			  (i + 1) * (2 * x_kol - 1) + j, (i + 1)* (2 * x_kol - 1) + j + 1, (i + 1)* (2 * x_kol - 1) + j + 2,
//			  (i + 2) * (2 * x_kol - 1) + j, (i + 2)* (2 * x_kol - 1) + j + 1, (i + 2)* (2 * x_kol - 1) + j + 2,
//			};
//
//			el_id++;
//		}
//	}
//
//	kol_S1nodes = 2 * (2 * y_kol - 1) + (2 * x_kol - 1);
//	S1.resize(kol_S1nodes);
//	int S1_id = 0;
//	// right
//	for (int i = 0; i < 2 * y_kol - 1; i++)
//	{
//		S1[S1_id] = i * (2 * x_kol - 1);
//		S1_id++;
//	}
//	// top
//	for (int i = 1; i < 2 * x_kol; i++)
//	{
//		S1[S1_id] = (2 * y_kol - 1) * (2 * x_kol - 1) - i;
//		S1_id++;
//	}
//	// left
//	for (int i = 1; i < 2 * y_kol; i++)
//	{
//		S1[S1_id] = i * (2 * x_kol - 1) - 1;
//		S1_id++;
//	}
//
//}
//int Mesh::MakeBiquadratic_test(string path)
//{
//	string file = path + "mesh.txt";
//	ifstream inf;
//	double tmpdbl;
//	int tmpint = 0, x_kol, y_kol;
//	double x_k, y_k, hx, hy, t1, t2;
//	vector<double> x, y;
//
//	inf.open(file);
//	inf >> t1 >> t2 >> x_kol >> x_k;
//
//	x.resize(2 * x_kol - 1);
//	int ii = 2;
//	x[0] = t1;
//
//	for (int i = 1; i < 2 * x_kol - 2; i += 2)
//	{
//		tmpdbl = x[i - 1] + hx;
//		x[i] = (tmpdbl + x[i - 1]) / 2.0;
//		x[i + 1] = tmpdbl;
//	}
//	inf.close();
//	inf.clear();
//
//	y.resize(2 * y_kol - 1);
//	y[0] = 0;
//	ii = 2;
//	for (int i = 1; i < 2 * y_kol - 2; i += 2, ii+=2)
//	{
//		tmpdbl = ii;
//		y[i] = (tmpdbl + y[i - 1]) / 2.0;
//		y[i + 1] = tmpdbl;
//	}
//	inf.close();
//	inf.clear();
//
//	int node_id = 0;
//	kol_nodes = (2 * y_kol - 1) * (2 * x_kol - 1);
//	nodes.resize(kol_nodes);
//
//	for (int i = 0; i < 2 * y_kol - 1; i++)
//	{
//		for (int j = 0; j < 2 * x_kol - 1; j++)
//		{
//			nodes[node_id].x = x[j];
//			nodes[node_id].y = y[i];
//			node_id++;
//		}
//	}
//
//	int el_id = 0;
//	for (int i = 0; i < 2 * y_kol - 2; i += 2)
//	{
//		for (int j = 0; j < 2 * x_kol - 2; j += 2)
//		{
//			elements[el_id].node_loc =
//			{ i * (2 * x_kol - 1) + j,      i * (2 * x_kol - 1) + j + 1,      i * (2 * x_kol - 1) + j + 2,
//			  (i + 1) * (2 * x_kol - 1) + j, (i + 1) * (2 * x_kol - 1) + j + 1, (i + 1) * (2 * x_kol - 1) + j + 2,
//			  (i + 2) * (2 * x_kol - 1) + j, (i + 2) * (2 * x_kol - 1) + j + 1, (i + 2) * (2 * x_kol - 1) + j + 2,
//			};
//
//			el_id++;
//		}
//	}
//	elements.resize(el_id);
//	kol_elements = el_id;
//
//	kol_S1nodes = 2 * (2 * y_kol - 1) + (2 * x_kol - 1);
//	S1.resize(kol_S1nodes);
//	int S1_id = 0;
//	// right
//	for (int i = 0; i < 2 * y_kol - 1; i++)
//	{
//		S1[S1_id] = i * (2 * x_kol - 1);
//		S1_id++;
//	}
//	// top
//	for (int i = 1; i < 2 * x_kol; i++)
//	{
//		S1[S1_id] = (2 * y_kol - 1) * (2 * x_kol - 1) - i;
//		S1_id++;
//	}
//	// left
//	for (int i = 1; i < 2 * y_kol; i++)
//	{
//		S1[S1_id] = i * (2 * x_kol - 1) - 1;
//		S1_id++;
//	}
//	return 0;
//}

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
	PointsOnAxis(in, all_R, count_r);
	PointsOnAxis(in, all_Z, count_z);
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

//int Mesh::MakeLinear_test(string path)
//{
//	string file = path + "r.dat";
//	ifstream inf;
//	double tmpdbl;
//	int tmpint = 0, x_kol, y_kol;
//	vector<double> x, y;
//
//	x_kol = 7;
//	y_kol = 7;
//
//	x.resize(x_kol);
//	for (int i = 0; i < x_kol; i++)
//	{
//		x[i] = 3.0 + i;
//	}
//	inf.close();
//	inf.clear();
//
//	y.resize(y_kol);
//	for (int i = 0; i < y_kol; i++)
//	{
//		y[i] = 2.0 + i;
//	}
//	inf.close();
//	inf.clear();
//
//	int node_id = 0;
//	kol_nodes = (y_kol) * (x_kol);
//	nodes.resize(kol_nodes);
//	for (int i = 0; i < y_kol; i++)
//	{
//		for (int j = 0; j < x_kol; j++)
//		{
//			nodes[node_id].x = x[j];
//			nodes[node_id].y = y[i];
//			node_id++;
//		}
//	}
//
//	int el_id = 0;
//	elements.resize((y_kol -1) * (x_kol - 1));
//
//	for (int i = 0; i < y_kol - 1; i++)
//	{
//		for (int j = 0; j < x_kol - 1; j++)
//		{
//			elements[el_id].node_loc =
//			{  i      * (x_kol) + j,      i  * (x_kol) + j + 1,
//			  (i + 1) * (x_kol) + j, (i + 1) * (x_kol) + j + 1
//			};
//			elements[el_id].material = 0;
//			el_id++;
//		}
//	}
//	kol_elements = el_id;
//
//	materials.resize(1);
//	materials[0].mu = 1;
//	materials[0].J = 0;
//
//
//	kol_S1nodes = 2 * (y_kol) + 2*(x_kol) + 8;
//	S1.resize(kol_S1nodes);
//	int S1_id = 8;
//	S1[0] = 16;
//	S1[1] = 17;
//	S1[2] = 18;
//	S1[3] = 23;
//	S1[4] = 25;
//	S1[5] = 30;
//	S1[6] = 31;
//	S1[7] = 32;
//
//
//	// right
//	for (int i = 0; i < y_kol; i++)
//	{
//		S1[S1_id] = i * (x_kol);
//		S1_id++;
//	}
//	// top
//	for (int i = 1; i < x_kol + 1; i++)
//	{
//		S1[S1_id] = (y_kol) * (x_kol) - i;
//		S1_id++;
//	}
//	// left
//	for (int i = 1; i < y_kol + 1; i++)
//	{
//		S1[S1_id] = i * (x_kol) - 1;
//		S1_id++;
//	}
//	// bottom
//	for (int i = 0; i < x_kol; i++)
//	{
//		S1[S1_id] = i;
//		S1_id++;
//	}
//	return 0;
//}


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