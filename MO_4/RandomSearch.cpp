#include "RandomSearch.h"

SimpleRandomSearch::SimpleRandomSearch() :
	gen(rd())
{
	x_a = 0;
	x_b = 0;
	y_a = 0;
	y_b = 0;
	eps = 1e-3;
	N_simple = 0;
	N = 0;
	m = 2;
	m_global = 1;
	func_0 = 0;
	x = { 0, 0 };
	a = { 0, 0 };
	b = { 0, 0 };
	a_start = { 0, 0 };
	b_start = { 0, 0 };
}

void SimpleRandomSearch::input()
{
	std::fstream in("in.txt");
	if (in.is_open())
	{
		in >> x_a;
		in >> x_b;
		in >> y_a;
		in >> y_b;
		in >> eps;
		in >> N_simple;
		in >> N;
		in >> m;
		in >> m_global;
		in >> a[0];
		in >> b[0];
		in >> a[1];
		in >> b[1];
		a_start[0] = a[0];
		b_start[0] = b[0];
		a_start[1] = a[1];
		b_start[1] = b[1];
	}
}

double SimpleRandomSearch::getRandomDouble(double min, double max)
{
	std::uniform_real_distribution<double> dist(min, max);
	return dist(gen);
}

void SimpleRandomSearch::Search(double x_min, double x_max, double y_min, double y_max, bool Ifbreak, bool IfStartPoint)
{
	double func_1 = 0;

	if (IfStartPoint)
	{
		x[0] = getRandomDouble(x_min, x_max);
		x[1] = getRandomDouble(y_min, y_max);

		//logPoint("S", x[0], x[1]);
	}
	func_0 = func.f(x[0], x[1]);

	for (int i = 0; i < N_simple; i++)
	{
		double x_1 = getRandomDouble(x_min, x_max);
		double y_1 = getRandomDouble(y_min, y_max);
		func_1 = func.f(x_1, y_1);

		//logAllPoint("P", x_1, y_1);

		//logPoint("R", x_1, y_1);

		if (func_1 < func_0)
		{
			x[0] = x_1;
			x[1] = y_1;
			func_0 = func_1;

			//logPoint("B", x[0], x[1]);

			if (Ifbreak) break;
		}
	}
}

void SimpleRandomSearch::Search()
{
	Search(x_a, x_b, y_a, y_b, false, true);
}

void SimpleRandomSearch::hyperquadrate(std::vector<double> a, std::vector<double> b, bool NeedStartPoint)
{
	double func_1 = 0;

	if (NeedStartPoint)
	{
		x[0] = getRandomDouble(a[0], b[0]);
		x[1] = getRandomDouble(a[1], b[1]);
	}
	func_0 = func.f(x[0], x[1]);

	for (int k = 0; k < N; k++)
	{
		logQuadrate("Q", a[0], b[0], a[1], b[1]);
		for (int i = 0; i < m; i++)
		{
			double x_1 = getRandomDouble(a[0], b[0]);
			double y_1 = getRandomDouble(a[1], b[1]);
			func_1 = func.f(x_1, y_1);

			logAllPoint("P", x_1, y_1);

			if (func_1 < func_0)
			{
				x[0] = x_1;
				x[1] = y_1;
				func_0 = func_1;

				logPoint("B", x[0], x[1]);
			}
		}

		double a_i = 0;
		double b_i = 0;

		for (int i = 0; i < 2; i++)
		{
			a_i = a[i];
			b_i = b[i];
			a[i] = x[i] - (b_i - a_i) / 2;
			b[i] = x[i] + (b_i - a_i) / 2;

			if (i == 0)
			{
				if (a[i] < x_a) a[i] = x_a;
				if (b[i] > x_b) b[i] = x_b;
			}
			else
			{
				if (a[i] < y_a) a[i] = y_a;
				if (b[i] > y_b) b[i] = y_b;
			}
		}
	}
	logPoint("L", x[0], x[1]);
}

void SimpleRandomSearch::hyperquadrate()
{
	hyperquadrate(a_start, b_start, true);
}

void SimpleRandomSearch::algorithm_1()
{
	double func_temp = 1e10;
	double x_temp = 0;
	double y_temp = 0;

	for (int i = 0; i < m_global; i++)
	{
		hyperquadrate(a_start, b_start, true);

		if (func_0 < func_temp)
		{
			func_temp = func_0;
			x_temp = x[0];
			y_temp = x[1];

			logPoint("B", x[0], x[1]);
		}
	}

	func_0 = func_temp;
	x[0] = x_temp;
	x[1] = y_temp;

	logPoint("L", x[0], x[1]);
}

void SimpleRandomSearch::algorithm_2()
{
	hyperquadrate(a_start, b_start, true);

	double best_f = func_0;
	std::vector<double> best_x = x;

	for (int k = 0; k < m_global; k++)
	{
		bool found = false;

		for (int attempt = 0; attempt < m; attempt++)
		{

			double x_rand = getRandomDouble(x_a, x_b);
			double y_rand = getRandomDouble(y_a, y_b);

			double f_rand = func.f(x_rand, y_rand);

			logPoint("R", x_rand, y_rand);

			if (f_rand < func_0)
			{
				x[0] = x_rand;
				x[1] = y_rand;
				func_0 = f_rand;

				logPoint("B", x[0], x[1]);

				found = true;
				break;
			}
		}

		if (!found)
			break;

		std::vector<double> a_local =
		{
			 std::max(x[0] - 1.0, x_a),
			 std::max(x[1] - 1.0, y_a)
		};
		std::vector<double> b_local =
		{
			 std::min(x[0] + 1.0, x_b),
			 std::min(x[1] + 1.0, y_b)
		};

		hyperquadrate(a_local, b_local, false);

		if (func_0 < best_f)
		{
			best_f = func_0;
			best_x = x;
		}
	}

	x = best_x;
	func_0 = best_f;
}

void SimpleRandomSearch::algorithm_3()
{
	x[0] = getRandomDouble(x_a, x_b);
	x[1] = getRandomDouble(y_a, y_b);

	logPoint("S", x[0], x[1]);

	hyperquadrate(a_start, b_start, false);

	double best_f = func_0;
	std::vector<double> best_x = x;

	for (int k = 0; k < m_global; k++)
	{
		logPoint("L", x[0], x[1]);

		double dx = getRandomDouble(-1.0, 1.0);
		double dy = getRandomDouble(-1.0, 1.0);

		double step = 5.0;

		std::vector<double> x_jump =
		{
			 x[0] + step * dx,
			 x[1] + step * dy
		};

		if (x_jump[0] < x_a || x_jump[0] > x_b ||
			x_jump[1] < y_a || x_jump[1] > y_b)
		{
			continue;
		}

		logPoint("J", x_jump[0], x_jump[1]);

		x = x_jump;
		func_0 = func.f(x[0], x[1]);

		std::vector<double> a_local =
		{
			 std::max(x[0] - 1.0, x_a),
			 std::max(x[1] - 1.0, y_a)
		};

		std::vector<double> b_local =
		{
			 std::min(x[0] + 1.0, x_b),
			 std::min(x[1] + 1.0, y_b)
		};

		hyperquadrate(a_local, b_local, false);

		if (func_0 < best_f)
		{
			best_f = func_0;
			best_x = x;

			logPoint("B", x[0], x[1]);
		}
		else
		{

			x = best_x;
			func_0 = best_f;
		}
	}

	x = best_x;
	func_0 = best_f;

	logPoint("F", x[0], x[1]);
}

void SimpleRandomSearch::logPoint(const std::string& type, double x, double y)
{
	if (log_file.is_open())
		log_file << type << " " << x << " " << y << "\n";
}

void SimpleRandomSearch::logQuadrate(const std::string& type, double x_min, double x_max, double y_min, double y_max)
{
	if (log_file.is_open())
	{
		log_file << type << " " << x_min << " " << x_max << " " << y_min << " " << y_max << "\n";
	}
}

void SimpleRandomSearch::logAllPoint(const std::string& type, double x, double y)
{
	if (log_file.is_open())
		log_file << type << " " << x << " " << y << "\n";
}