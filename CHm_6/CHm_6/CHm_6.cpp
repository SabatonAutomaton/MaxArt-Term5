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
int* fict = nullptr;
int* kray = nullptr;
int n_x = 0, n_y = 0, m = 0, maxiter = 10, countDiag = 5, maindiag = 0;
double w = 1;//параметр w
double accuracy = 1e-5;
double a_x = 0, b_x = 0, h_x = 0, q_x = 0;
double a_y = 0, b_y = 0, h_y = 0, q_y = 0;

double lambda = 1;
double gamma = 1;
double thetha = 1;

vector<double> mesh_x;
vector<double> mesh_y;

double func_kray(double x, double y, int choise)
{
	switch (choise)
	{
	case 1:
		return x+y;
	case 2:
		return 1;
	default:
		return 0;
	}
}


static void read(int choise)
{
	ifstream f("mesh.txt");
	ifstream fp("params.txt");
	ifstream ff("fict.txt");
	ifstream fk("kray.txt");

	if (!f || !fp || !ff || !fk)
	{
		cerr << "File open error" << endl;
		return;
	}

	fp >> accuracy;
	fp >> maxiter;
	fp >> w;

	f >> n_x;//количество точек по х
	f >> q_x;
	f >> a_x;
	f >> b_x;



	f >> n_y;//количество точек по y
	f >> q_y;
	f >> a_y;
	f >> b_y;

	double hi = 0;
	double sum = 0;

	switch (choise)
	{
	case 0:
		h_x = (b_x - a_x) / (n_x - 1);
		h_y = (b_y - a_y) / (n_y - 1);
		for (int i = 0; i < n_x; i++)
			mesh_x.push_back(a_x + i * h_x);
		for (int i = 0; i < n_y; i++)
			mesh_y.push_back(a_y + i * h_y);
		break;
	case 1:
		h_x = (b_x - a_x) * (1 - q_x) / (1 - pow(q_x, n_x - 1));
		h_y = (b_y - a_y) * (1 - q_y) / (1 - pow(q_y, n_y - 1));

		hi = h_x;
		mesh_x.push_back(a_x);
		for (int i = 1; i < n_x - 1; i++)
		{
			sum += h_x * pow(q_x, i - 1);
			mesh_x.push_back(a_x + sum);

		}
		mesh_x.push_back(b_x);

		hi = h_y;
		sum = 0;
		mesh_y.push_back(a_y);
		for (int i = 1; i < n_y - 1; i++)
		{
			sum += h_y * pow(q_y, i - 1);
			mesh_y.push_back(a_y + sum);
		}
		mesh_y.push_back(b_y);

		break;
	default:
		break;
	}




	initialApprox = new double[n_x * n_y] {0};
	new_initialApprox = new double[n_x * n_y]();
	result = new int[n_x * n_y]();

	maindiag = 2;

	fict = new int[n_x * n_y] {0};
	int x_a = 0, x_b = 0, y_a = 0, y_b = 0, type = 0;

	while (ff >> type >> x_a >> x_b >> y_a >> y_b)
	{
		for (int i = x_a; i <= x_b; i++)
			for (int j = y_a; j <= y_b; j++)
			{
				if (!type)
					fict[j * n_x + i] = 1;
				else if (fict[j * n_x + i] == 1)
					fict[j * n_x + i] = 2;
				else if (fict[j * n_x + i] != 2)
					fict[j * n_x + i] = 0;
			}
	}

	kray = new int[n_x * n_y] {0};
	x_a = 0, x_b = 0, y_a = 0, y_b = 0, type = 0;
	int number_func = 0;

	while (fk >> type >> number_func >> x_a >> x_b >> y_a >> y_b)
	{

		for (int i = x_a; i <= x_b; i++)
		{
			for (int j = y_a; j <= y_b; j++)
			{
				int k = n_x * j + i;
				if (type == 1)
					kray[k] = 1;
				else if (type == 2)
					if (kray[k]!=1)
						kray[k] = number_func;
			}
		}
	}
	//отладочная
	for (int i = 0; i < n_x * n_y; i++)
		cout << "kray[" << i << "] = " << kray[i] << " ";
	cout << '\n';
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

	vect = new double[n_x * n_y];

	for (int i = 0; i < n_x; i++)
	{
		if (kray[i] == 1)
		{
			A[2][i] = 1;
			vect[i] = func_kray(mesh_x[i], mesh_y[0], kray[i]);
		}
		else
		{
			double h = mesh_y[1] - mesh_y[0];
			A[2][i] = -lambda / h;
			A[4][i] = lambda / h;
			vect[i] = func_kray(mesh_x[i], mesh_y[0], kray[i]);
		}

		if (kray[n_x * n_y - i - 1] == 1)
		{
			A[2][n_x * n_y - i - 1] = 1;
			vect[n_x * n_y - i - 1] = func_kray(mesh_x[n_x - i - 1], mesh_y[n_y - 1], kray[n_x * n_y - i - 1]);
		}
		else
		{
			double h = mesh_y[n_y - 1] - mesh_y[n_y - 2];
			A[0][n_x * n_y - i - 1 + ind[0]] = -lambda / h;
			A[2][n_x * n_y - i - 1] = lambda / h;
			vect[n_x * n_y - i - 1] = func_kray(mesh_x[n_x - i - 1], mesh_y[n_y - 1], kray[n_x * n_y - i - 1]);
		}

	}

	for (int i = 1; i < n_y - 1; i++)
	{
		double hn_y = mesh_y[i + 1] - mesh_y[i - 1];

		if (kray[n_x * i] == 1)
		{
			A[2][n_x * i] = 1;
			vect[n_x * i] = func_kray(mesh_x[0], mesh_y[i], kray[n_x * i]);
		}
		else
		{
			double h = mesh_x[1] - mesh_x[0];
			A[2][n_x * i] = -lambda / h;
			A[3][n_x * i] = lambda / h;
			vect[n_x * i] = func_kray(mesh_x[0], mesh_y[i], kray[n_x * i]);
		}

		if (kray[n_x * i + n_x - 1] == 1)
		{
			A[2][n_x * i + n_x - 1] = 1;
			vect[n_x * i + n_x - 1] = func_kray(mesh_x[n_x - 1], mesh_y[i], kray[n_x * i + n_x - 1]);
		}
		else
		{
			double h = mesh_x[n_x - 1] - mesh_x[n_x - 2];
			A[1][n_x * i + n_x - 1 + ind[1]] = -lambda / h;
			A[2][n_x * i + n_x - 1] = lambda / h;
			vect[n_x * i + n_x - 1] = func_kray(mesh_x[n_x - 1], mesh_y[i], kray[n_x * i + n_x - 1]);
		}

		for (int j = 1; j < n_x - 1; j++)
		{
			int k = n_x * i + j;
			double hn_x = mesh_x[j + 1] - mesh_x[j - 1];

			if (fict[k] == 2)
			{
				if (kray[k] == 1)
					A[2][k + ind[2]] = 1;
				vect[k] = func_kray(mesh_x[j], mesh_y[i], kray[k]);
			}
			else if (fict[k] == 0)
			{
				A[0][k + ind[0]] = -lambda / (hn_y * hn_y);
				A[1][k + ind[1]] = -lambda / (hn_x * hn_x);
				A[2][k + ind[2]] = 2 * lambda * (1 / (hn_x * hn_x) + 1 / (hn_y * hn_y)) + gamma;
				A[3][k] = -lambda / (hn_x * hn_x);
				A[4][k] = -lambda / (hn_y * hn_y);

				vect[k] = func_kray(mesh_x[j], mesh_y[i], 1);
			}
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
	{
		if (fict[i] != 1)
			out << fixed << setprecision(16) << initialApprox[i] << '\n';
		else
			out << fixed << setprecision(16) << 0 << '\n';
	}
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
	bool choise = 0;
	cout << "Enter type mesh(0,1) 0-Uniform, 1-Not uniform\n";
	cin >> choise;

	read(choise);
	int n = n_x * n_y;

	build_portrait();
	build_matrix();

	out_console(n);


	//cout << "Enter method(0,1) 0-Jacobi, 1-Gauss_Seidel\n";
	//cin >> choise;
	choise = 1;
	cout << "k     residual       w\n";
	Jacobi_Gauss_Seidel(n, choise, w);


	out_file(n);


	for (int i = 0; i < countDiag; i++)
	{
		delete[] A[i];
	}
	delete[] A;
	delete[] ind;
	delete[] vect;
	delete[] fict;
}