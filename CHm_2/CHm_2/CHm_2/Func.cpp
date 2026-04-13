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

	case 7:
		return x * x;

	case 8:
		return x * x;

	case 9:
		return x + currentTime;

	case 10:
		return x * x + currentTime;

	case 11:
		return x + currentTime * currentTime;

	case 12:
		return x + currentTime;

	case 13:
		return x * x + currentTime * currentTime;

	case 14:
		return x + currentTime;

	case 15:
		return x * x + currentTime;

	case 16:
		return x + currentTime;

	case 17:
		return x * x + currentTime * currentTime;

	case 18:
		return x * x * x;

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

	case 7:
		return -12.0 * x * x - 10.0;

	case 8:
		return 0.8 * x * x * sin( x * x / 5.0 ) - 2.0 * cos( x * x / 5.0 );

	case 9:
		return 1.0;

	case 10:
		return -1.0;

	case 11:
		return 2.0 * currentTime;

	case 12:
		return 0.0;

	case 13:
		return 2.0 * currentTime - 6.0 * x * x - 2.0 * currentTime * currentTime;

	case 14:
		return -1.0;

	case 15:
		return -12.0 * x * x - 4.0 * currentTime - 9.0;

	case 16:
		return 1.0 + 0.2 * sin( ( x + currentTime ) / 5.0 );

	case 17:
		return 2.0 * currentTime
			- 2.0 * cos( ( x * x + currentTime * currentTime ) / 5.0 )
			+ 0.8 * x * x * sin( ( x * x + currentTime * currentTime ) / 5.0 );

	case 18:
		return -6.0 * x;

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
		return u;
	case 5:
		return u;

	case 6:
		return cos( u );

	case 7:
		return 2.0 * u + 5.0;

	case 8:
		return cos( u / 5.0 );

	case 9:
	case 10:
	case 11:
		return 1.0;

	case 12:
	case 13:
		return u;

	case 14:
		return 2.0 * u + 5.0;

	case 15:
		return 2.0 * u + 5.0;

	case 16:
	case 17:
		return cos( u / 5.0 );

	case 18:
		return 1.0;

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
		return 1.0;
	case 5:
		return 1.0;

	case 6:
		return -sin( u );

	case 7:
		return 2.0;

	case 8:
		return -sin( u / 5.0 ) / 5.0;

	case 9:
	case 10:
	case 11:
		return 0.0;

	case 12:
	case 13:
		return 1.0;

	case 14:
		return 2.0;

	case 15:
		return 2.0;

	case 16:
	case 17:
		return -sin( u / 5.0 ) / 5.0;

	case 18:
		return 0.0;

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
