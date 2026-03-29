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
		return x * x * x;

	case 4:
		return x * x * x * x;

	case 5:
		return x * x * x * x * x;

	case 6:
		return x / (1.0 + x);

	case 7:
		return sin(std::numbers::pi * x);

	case 8:
		return exp(-currentTime) * sin(std::numbers::pi * x);

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
		return -6.0 * x;

	case 4:
		return -12.0 * x * x;

	case 5:
		return -20.0 * x * x * x;

	case 6:
		return 2.0 / pow(1.0 + x, 3);

	case 7:
		return std::numbers::pi * std::numbers::pi * sin(std::numbers::pi * x);

	case 8:
		return (std::numbers::pi * std::numbers::pi - 1.0) * exp(-currentTime) * sin(std::numbers::pi * x);

	default:
		return 0.0;
	}
}

double Functions::lambda(double u)
{
	switch (id)
	{
	default:
		return 1.0;
	}
}

double Functions::dlambda(double u)
{
	switch (id)
	{
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
