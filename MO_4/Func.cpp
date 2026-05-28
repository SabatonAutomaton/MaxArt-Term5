#include "Func.h"

double Function::f(double x, double y)
{
	count_use_f++;
	double sum = 0;

	for (int i = 0; i < 6; i++)
	{
		sum += -C[i] / (1 + (x - a[i]) * (x - a[i]) + (y - b[i]) * (y - b[i]));
	}
	return sum;
}
