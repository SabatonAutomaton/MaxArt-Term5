#pragma once

#include "Func.h"
#include <fstream>
#include <random>
#include <vector>
#include <string>

class SimpleRandomSearch
{
private:

	double x_a;
	double x_b;
	double y_a;
	double y_b;
	double eps;
	int N_simple;
	int N;
	int m;
	int m_global;
	
	std::random_device rd;
	std::mt19937 gen;

	std::vector<double> a;
	std::vector<double> b;
	std::vector<double> a_start;
	std::vector<double> b_start;

public:
	double func_0;
	std::vector<double> x;
	Function func;
	std::ofstream log_file;

public:
	void input();
	double getRandomDouble(double min, double max);
	void Search(double x_min, double x_max, double y_min, double y_max, bool Ifbreak, bool IfStartPoint);
	void Search();
	SimpleRandomSearch();
	void hyperquadrate();
	void hyperquadrate(std::vector<double> a, std::vector<double> b, bool NeedStartPoint);
	void algorithm_1();
	void algorithm_2();
	void algorithm_3();
	void logPoint(const std::string& type, double x, double y);
	void logQuadrate(const std::string& type, double x_min, double x_max, double y_min, double y_max);
	void logAllPoint(const std::string& type, double x, double y);
};

