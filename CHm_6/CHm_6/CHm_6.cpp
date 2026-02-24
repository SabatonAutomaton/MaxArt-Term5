#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

double** A = nullptr;
int* ind = nullptr;
double* vect = nullptr;
double* initialApprox = nullptr;
double* new_initialApprox = nullptr;
int* result = nullptr;
int n_x = 0, n_y = 0, m = 0, maxiter = 10, countDiag = 5, maindiag = 0;
double w = 1;//параметр w
double accuracy = 1e-5;
double a_x = 0, b_x = 0, h_x = 0;
double a_y = 0, b_y = 0, h_y = 0;

double lambda = 1;
double gamma = 1;

vector<double> mesh_x;
vector<double> mesh_y;

double u(double x, double y)
{
	return 1;
}

static void read()
{
	ifstream f("mesh.txt");
	ifstream fv("vector.txt");
	ifstream fp("params.txt");

	if (!f || !fv || !fp)
	{
		cerr << "File open error" << endl;
		return;
	}

	fp >> accuracy;
	fp >> maxiter;
	fp >> w;

	f >> n_x;//количество точек по х
	f >> a_x;
	f >> b_x;

	h_x = (b_x - a_x) / (n_x - 1);

	f >> n_y;//количество точек по х
	f >> a_y;
	f >> b_y;

	h_y = (b_y - a_y) / (n_y - 1);

	//vector right part
	vect = new double[n_x * n_y];
	for (int i = 0; i < n_x * n_y; i++)
		fv >> vect[i];

	for (int i = 0; i < n_x; i++)
		mesh_x.push_back(a_x + i * h_x);
	for (int i = 0; i < n_y; i++)
		mesh_y.push_back(a_y + i * h_y);

	initialApprox = new double[n_x * n_y] {0};
	new_initialApprox = new double[n_x * n_y]();
	result = new int[n_x * n_y]();

	maindiag = 2;
}

void build_portrait()
{
	countDiag = 5;
	A = new double* [countDiag];
	for (int i = 0; i < countDiag; i++)
		A[i] = new double[n_x * n_y] {0};

	m = n_x - 2;
	ind = new int[countDiag] {-m - 2, -1, 0, 1, 2 + m};
}

void build_matrix()
{
	for (int i = 0; i < n_x; i++)
	{
		A[2][i] = 1;
		A[2][n_x * n_y - i - 1] = 1;

		vect[i] = u(mesh_x[i], mesh_y[0]);
		vect[n_x * n_y - i - 1] = u(mesh_x[n_x - i - 1], mesh_y[n_y - 1]);
	}

	for (int i = 1; i < n_y - 1; i++)
	{
		A[2][n_x * i] = 1;
		A[2][n_x * i + n_x - 1] = 1;

		vect[n_x * i] = u(mesh_x[0], mesh_y[i]);
		vect[n_x * i + n_x - 1] = u(mesh_x[n_x - 1], mesh_y[i]);

		for (int j = 1; j < n_x - 1; j++)
		{
			int k = n_x * i + j;
			A[0][k + ind[0]] = -lambda / (h_y * h_y);
			A[1][k + ind[1]] = -lambda / (h_x * h_x);
			A[2][k + ind[2]] = 2 * lambda * (1 / (h_x * h_x) + 1 / (h_y * h_y)) + gamma;
			A[3][k] = -lambda / (h_x * h_x);
			A[4][k] = -lambda / (h_y * h_y);
		}
	}



}

void out_console(int n) {

	for (int i = 0; i < countDiag; i++)
		cout << "i[" << i << "] = " << ind[i] << " ";
	cout << '\n';

	for (int i = 0; i < countDiag; i++)
	{
		for (int j = 0; j < n; j++)
			cout << "A[" << i << "][" << j << "] = " << A[i][j] << " ";
		cout << '\n';
	}

	for (int i = 0; i < n; i++)
		cout << "vect[" << i << "] = " << vect[i] << " ";
	cout << '\n';

	for (int i = 0; i < n; i++)
		cout << "x[" << i << "] = " << fixed << setprecision(5) << initialApprox[i] << " ";
	cout << '\n';

}

void out_file(int n)
{
	fstream out("out.txt", 2);
	for (int i = 0; i < n; i++)
		out << fixed << setprecision(16) << initialApprox[i] << '\n';
}

void multiplice(int n)
{
	for (int i = 0; i < n; i++)
		new_initialApprox[i] = 0;

	for (int i = 0; i < n; i++)
	{
		for (int c = 0; c < countDiag; c++)// бегаем "c" по ind[c] если диагональ выше строки или ниже но её элемент на нашей строке , то умножаем
		{
			int j = i + ind[c];
			if (j >= 0 && j < n)
			{
				if (ind[c] >= 0)// выше строки
				{
					new_initialApprox[i] += A[c][i] * initialApprox[j];
				}
				else// совпадает элемент диагонали (ниже строки)
				{
					new_initialApprox[i] += A[c][j] * initialApprox[j];
				}
			}
		}
	}
}

double residual(int n)
{
	multiplice(n);
	double norm_1 = 0, norm_2 = 0;
	for (int i = 0; i < n; i++)
	{
		new_initialApprox[i] = vect[i] - new_initialApprox[i];
		norm_1 += new_initialApprox[i] * new_initialApprox[i];
		norm_2 += vect[i] * vect[i];
	}
	norm_1 = sqrt(norm_1);
	norm_2 = sqrt(norm_2);

	return norm_1 / norm_2;
}

void Jacobi_Gauss_Seidel(int n, bool choise, double w)
{
	double* test = nullptr;
	if (!choise)
		test = initialApprox;
	else
		test = new_initialApprox;

	for (int k = 1; k <= maxiter; k++)// Максимальное количество итераций
	{
		for (int i = 0; i < n; i++)
		{
			double sum1 = 0, sum2 = 0;// счёт суммы

			for (int c = 0; c < countDiag; c++)// бегаем "c" по ind[c] если диагональ выше строки или ниже но её элемент на нашей строке , то умножаем
			{
				int j = i + ind[c];
				if (j >= 0 && j < n)
				{
					if (ind[c] >= 0)// выше строки или на ней
					{
						sum1 += A[c][i] * initialApprox[j];
					}
					else// совпадает элемент диагонали (ниже строки)
					{
						sum2 += A[c][j] * test[j];
					}
				}
			}
			new_initialApprox[i] = initialApprox[i] + (w / A[maindiag][i]) * (vect[i] - sum1 - sum2);
		}
		for (int i = 0; i < n; i++)
			initialApprox[i] = new_initialApprox[i];

		// Добавлен выход по относительной невязке

		double res = residual(n);
		cout << "k = " << k << "; residual = " << fixed << setprecision(16) << res << defaultfloat << '\n';

		if (res < accuracy)
		{
			cout << "residual = " << fixed << setprecision(16) << res << " < " << accuracy << '\n' << defaultfloat;
			return;
		}
	}
}

int main()
{

	read();
	int n = n_x * n_y;

	build_portrait();
	build_matrix();

	out_console(n);

	bool choise = 0;
	cout << "Enter method(0,1) 0-Jacobi, 1-Gauss_Seidel\n";
	cin >> choise;

	cout << "k     residual       w\n";
	Jacobi_Gauss_Seidel(n,choise,w);


	out_file(n);


	for (int i = 0; i < countDiag; i++)
	{
		delete[] A[i];
	}
	delete[] A;
	delete[] ind;
	delete[] vect;
}
