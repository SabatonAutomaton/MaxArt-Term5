#include "Func.h"
#include <numbers>

double Functions::u(double x)
{
	switch (id)
	{
	case 1:
		return x;

	case 2:
		return x * x;

	case 3:
		return currentTime * currentTime;

	case 4:
		return x * x;

	case 5:
		return x * x + currentTime;

	case 6:
		return x * x;

	default:
		return x;
	}
}

double Functions::f(double x)
{
	switch (id)
	{
	case 1:
		return 0.0;

	case 2:
		return -2.0;

	case 3:
		return 2.0 * currentTime;

	case 4:
		return -6.0 * x * x;

	case 5:
		return 1.0 - 6.0 * x * x - 2.0 * currentTime;

	case 6:
		return 4.0 * x * x * sin( x * x ) - 2.0 * cos( x * x );

	default:
		return 0.0;
	}
}

double Functions::lambda(double u)
{
	switch (id)
	{
	case 1:
	case 2:
	case 3:
		return 1.0;

	case 4:
	case 5:
		return u;

	case 6:
		return cos( u );

	default:
		return 1.0;
	}
}

double Functions::dlambda(double u)
{
	switch (id)
	{
	case 1:
	case 2:
	case 3:
		return 0.0;

	case 4:
	case 5:
		return 1.0;

	case 6:
		return -sin( u );

	default:
		return 0.0;
	}
}

double Functions::sigma(double u)
{
	switch (id)
	{
	default:
		return 1.0;
	}
}

double Functions::sigmaTime(double t)
{
	switch (id)
	{
	default:
		return 1.0;
	}
}

double Functions::theta(double u)
{
	switch (id)
	{
	default:
		return 0.0;
	}
}

double Functions::beta(double u)
{
	switch (id)
	{
	default:
		return 0.0;
	}
}

double Functions::dbeta(double u)
{
	switch (id)
	{
	default:
		return 0.0;
	}
}

double Functions::uBeta(double u)
{
	switch (id)
	{
	default:
		return 0.0;
	}
}

double Functions::boundaryValue(int number, double x)
{
	switch (number)
	{
	case 1:
		return u(x);
	default:
		return 0.0;
	}
}
