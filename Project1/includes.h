#pragma once
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <set>
#include <sstream>


using namespace std; // отмена необходимости писать "std::" при работе с потоками

const int VARIABLE_COUNT = 4;
const double BIG_VALUE = 1e15;
const double MU0 = 4 * M_PI * 1e-7;
const int MAX_ITER = 100;
const double EPS_NONLINEAR = 1e-6;

template <typename T>
string NumberToString(T Number)
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}

typedef // создаем новый прототип (в данном случае указатель на функцию)
	double // возвращаемое значение (такое же как в функциях)
		(*basis_func) // имя прототипа (в коде употребляется без звездочки)
			(double); // список параметров (такое же как в функциях)

//---------------------------линейные одномерные-------------------------------------------
inline double linear1(double ksi)
{
	return 1 - ksi;
}
inline double linear2(double ksi)
{
	return ksi;
}
//--------------------------------------------------------------------------------


//---------------одномерные квадратичные функции, где ksi = (x - x_(2k-1)) / (x_(2k+1) - x_(2k-1)))-------------------------------------------------------
inline double quadratic1(double ksi)
{
	return 2 * (ksi - 0.5) * (ksi - 1);
}
inline double quadratic2(double ksi)
{
	return -4 * ksi * (ksi - 1);
}
inline double quadratic3(double ksi)
{
	return 2 * (ksi - 0.5) * ksi;
}
inline double getKsi(double x, double x1, double x2)
{
	return(x - x1) / (x2 - x1);
}
//--------------------------------------------------------------------------------

//--------------------------test--------------------------------------------------
inline double f_test(double x, double y)
{
	return 0;
}

inline double u_test(double x, double y)
{
	return x + y;
}

inline double mu_test()
{
	return 1;
}
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
const vector<vector<double>> G_linear = { 
	{1, -1}, 
	{-1, 1} }; // * 1/h
const vector<vector<double>> M_linear = { 
	{1.0 / 3.0, 1.0 / 6.0},
	{1.0 / 6.0, 1.0 / 3.0} }; // * h
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
const vector<vector<double>> G_quadratic = { 
	{ 7.0 / 3.0, -8.0 / 3.0,  1.0 / 3.0},
	{-8.0 / 3.0, 16.0 / 3.0, -8.0 / 3.0},
    { 1.0 / 3.0, -8.0 / 3.0,  7.0 / 3.0} }; // * 1/h
const vector<vector<double>> M_quadratic = {
	{ 4.0 / 30.0,  2.0 / 30.0, -1.0 / 30.0},
	{ 2.0 / 30.0, 16.0 / 30.0,  2.0 / 30.0},
	{-1.0 / 30.0,  2.0 / 30.0,  4.0 / 30.0} }; // * h

const vector<vector<double>> G_quadratic_r_1 = { 
	{7., -8., 1.},
	{-8., 16., -8.},
	{1., -8., 7.}};                            // G = 1/(3h)*(r_cur*G_quadratic_r_1 + h/2 * G_quadratic_r_2)
const vector<vector<double>> G_quadratic_r_2 = {
	{3., -4., 1.},
	{-4., 16., -12},
	{1., -12., 11} };                          // G = 1/(3h)*(r_cur*G_quadratic_r_1 + h/2 * G_quadratic_r_2)
const vector<vector<double>> M_quadratic_r_1 = {
	{4., 2., -1.},
	{2., 16., 2.},
	{-1., 2., 4.} };                           // M = h / 30 * (r_cur*M_quadratic_r_1 + h/2 * M_quadratic_r_2)
const vector<vector<double>> M_quadratic_r_2 = {
	{1., 0., -1.},
	{0., 16., 4.},
	{-1., 4., 7.} };                           // M = h / 30 * (r_cur*M_quadratic_r_1 + h/2 * M_quadratic_r_2)

inline int mu(int i)
{
	return (i % 3);
}
inline int nu(int i)
{
	return (i / 3);
}
//--------------------------------------------------------------------------------
