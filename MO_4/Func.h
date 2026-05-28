#pragma once
#include <vector>

class Function
{
private:
	std::vector<double> C = { 4, 9, 1, 7, 5, 6 };
	std::vector<double> a = { 7, -9, 6, -8, -10, -2 };
	std::vector<double> b = { 9, -1, 5, -2, -8, -4 };

public:
	double f(double x, double y);
	int count_use_f = 0;
};