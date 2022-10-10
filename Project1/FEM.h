#pragma once
#include "includes.h"
#include "cgm.h"
#include "Solver.h"
#include "Table.h"

struct Receiver
{
    double r[2];
    double z[2];
    double V, V_true;

    Receiver() { r[0] = 100; r[1] = 200; z[0] = 0; z[1] = 0; V = 0; V_true = 0; }
};

struct VEL
{
    double r[2];
    double z[2];
    double J[2];

    VEL() { r[0] = 0; z[0] = 0; r[1] = 0; z[1] = 1; J[0] = 1, J[1] = -1; }
};

struct Node
{
    double r;
    double z;

    Node() { r = 0; z = 0; }
    Node(double _r, double _z) { r = _r; z = _z; }
};

struct Elem
{
    std::vector<int> node_loc; // ���������� ������ ��������� ����� �����
    int material; // ����� ���������

    Elem()
    {
        node_loc.resize(9);
        material = 0;
    }
};

class Mesh
{
    int CreateMesh(string path);
    int ReadMat(string path);
public:
	vector<Elem> elements; // ��� �������� �����  sigma[elements[i].mat]
    vector<Node> nodes;         // ��� ���� � ������� ���������� ���������
    vector<double> sigma; // ��� ��������� �� �������� sigma[3]
    vector<int> S1;      // S1[j] �� j-�� ���� ������ ������� ������� 1 ����

    int kol_nodes, kol_elements, kol_S1nodes, kol_mat;

    int ReadMesh(string path);
    //int ReadMesh(string path);
    //int MakeBiquadratic(string path);
    //int MakeBiquadratic_test(string path);
    //int MakeLinear_test(string path);

    Node GetCentre(int el_id);

};

class SLAU
{

public:
    vector<int> ia, ja;
    //vector<double> di, a, b;
    vector<double> di, au, al, b;

    int GeneratePortret(Mesh& mesh);
    int AddLocal(vector<vector<double>>& A_loc, vector<double>& b_loc, int el_id, Mesh& mesh);
    int SetS1(Mesh& mesh);
    int SetS1_null(Mesh& mesh);

    //int SetS1_test(Mesh& mesh); // ���������� ������ �� �� �����, � �������� ������� �� includes 

    int ShowMatrix();   // ����������� ������� ������, ������� � ��������
    void ClearValues(); // ��������� 0 di, a, b

    int CheakG(Mesh& mesh);
};

//class FEM
//{
//    vector<vector<double>> A_loc;
//    vector<double> b_loc;
//    vector<double> q, q_prev;
//
//    int IsLinear = 1;
//
//    vector<int> conner_loc; // ��������� ������ ������� ����� ���������; 
//                            // [0] = ����� ������ [1] = ������ ������
//                            // [2] = ����� ������� [3] = ������ �������
//
//    vector<basis_func> basis; // �������� ���������� �������; ksi = (x - x_min) / (x_max - x_min))
//
//    int Init(string path);
//    int Init_biquadratic(string path);
//
//    int GetG_Loc(int el_id, double B = 0);
//    int GetG_Loc_test(int el_id);
//    int GetG_Loc_biquadratic(int el_id, double B = 0);
//    int GetG_Loc_biquadratic_test(int el_id);
//
//    int Getb_Loc(int el_id);
//    int Getb_Loc_test(int el_id);
//    int Getb_Loc_biquadratic(int el_id);
//    int Getb_Loc_biquadratic_test(int el_id);
//
//    double CheckRes();
//
//public:
//	Mesh mesh;
//    SLAU slau;
//    double w = 0.1; // �������� ���������� ��� ���������� ������
//
//    int SolveTask(string path);                  // ���������� �����
//    int SolveTask_test(string path);             // ���������� ������ �� �� �����, � �������� ������� �� includes 
//    int SolveTask_biquadratic(string path);      // �������������� �����
//    int SolveTask_biquadratic_test(string path); // ���������� ������ �� �� �����, � �������� ������� �� includes   
//
//    int GetB(vector<double> &Bx, vector<double> &By); // �������� (��, ��) �� ���� �����; ����������� ������� ��������������� ���� �����
//    int GetAz(vector<double> &Az);                    // �������� �z �� ���� �����
//
//    double GetAzInPoint(double x, double y, int cur_el = -1); // �������� �z � ����� (����� ������� �������, ���� �� ��������)
//    pair<double,double> GetBInPoint(double x, double y, double& Az, int cur_el = -1, double hx = 1e-3, double hy = 1e-3); // �������� Bx, By, Az � ����� (����� ������� �������, ���� �� ��������) (����� ������� ��� ��� �����������)
//};

class FEM_electro
{
    vector<vector<double>> A_loc;
    vector<double> b_loc;
    vector<double> q;
    int kol_rec, kol_vel;
    map<int, double> f_ist;  // [���� ����� ����] = ��������
    //vector<int> conner_loc; // ��������� ������ ������� ����� ���������; 
                            // [0] = ����� ������ [1] = ������ ������
                            // [2] = ����� ������� [3] = ������ �������

    vector<basis_func> basis; // �������� ���������� �������; ksi = (x - x_min) / (x_max - x_min))

    //int Init_biquadratic(string path);

    int GetG_Loc_biquadratic(int el_id);
    //int GetG_Loc_biquadratic_test(int el_id);

    int Getb_Loc_biquadratic(int el_id);
    //int Getb_Loc_biquadratic_test(int el_id);

    double GetF(double r, double z);
    double GetF(int node_id);


    int ReadData(string path);
    int FindElem(double r, double z);
public:

    Mesh mesh;
    SLAU slau;
    vector<Receiver> receivers;
    vector<VEL> VELs;

    int InitTask(string path);
    int DirectTask();
    int DirectTask(int mat_id, double h_sigma);
    double V_in_point(double r, double z);
    int V_in_rec(vector<double> &V);
    double GetUTest(double r, double z);

    int PrintField(string filename);
    int PrintLine(string filename, double z, double r); // V ����� ����� z, ������� � r

    double GetInverseFunc();

};