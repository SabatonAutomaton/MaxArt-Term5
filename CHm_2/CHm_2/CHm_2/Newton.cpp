#include "Newton.h"
#include <iostream>
#include <fstream>
#include <cmath>

double NewtonMethod::u(double x)
{
   return x * x;  // начальное приближение
}

double NewtonMethod::f(double x)
{
   return x * x - 2;
}

double NewtonMethod::Condition(int number, double x)
{
   switch (number)
   {
   case 1:
      return u(x);
      break;
   default:
      return 0;
   }
}

void NewtonMethod::ApplyBoundaryConditions()
{
   // Левая граница
   if (elements[0].cond == 1)
   {
      matrix.di[0] = 1.0;
      b[0] = Condition(1, 0.0);
      // Обнуляем связи с другими узлами
      matrix.ggu[0][0] = 0.0;
      matrix.ggu[0][1] = 0.0;
   }

   // Правая граница
   if (elements[elemCount - 1].cond == 1)
   {
      matrix.di[n - 1] = 1.0;
      b[n - 1] = Condition(1, 3.0);
      // Обнуляем связи с другими узлами
      matrix.ggl[n - 1][0] = 0.0;
      matrix.ggl[n - 1][1] = 0.0;
   }
}

void NewtonMethod::Input()
{
   mesh.InputMeshes();

   std::ifstream file("params.txt");
   if (!file)
   {
      std::cerr << "File open error" << std::endl;
      return;
   }
   file >> eps >> maxIter >> boundaryConditions[0] >> boundaryConditions[1];

   elemCount = mesh.n_x - 1;
   n = 3 * elemCount - elemCount + 1;  // для квадратичных элементов

   // Инициализация матрицы (ленточная)
   matrix.n = n;
   matrix.ggl.resize(n);
   matrix.ggu.resize(n);
   for (int i = 0; i < n; i++)
   {
      matrix.ggl[i] = std::vector<double>(2, 0.0);
      matrix.ggu[i] = std::vector<double>(2, 0.0);
   }
   matrix.di = std::vector<double>(n, 0.0);

   b.resize(n);

   // Инициализация начального приближения q
   q = std::vector<double>(n, 0.0);
   for (int i = 0; i <= elemCount; i++)
   {
      q[i * 2] = u(mesh.meshX[i]);
      if (i < elemCount)
         q[i * 2 + 1] = u((mesh.meshX[i] + mesh.meshX[i + 1]) / 2.0);
   }

   // Создание элементов
   std::vector<double> coords;
   elements.resize(elemCount);
   std::vector<double> localB = std::vector<double>(3, 0.0);

   for (int i = 0; i < elemCount; i++)
   {
      coords = { mesh.meshX[i], mesh.meshX[i + 1] };

      // Локальный вектор f (значения в узлах элемента)
      std::vector<double> localF(3);
      localF[0] = f(mesh.meshX[i]);
      localF[1] = f((mesh.meshX[i] + mesh.meshX[i + 1]) / 2.0);
      localF[2] = f(mesh.meshX[i + 1]);

      elements[i] = Element(coords, 0, q, localF);
   }

   elements[0].cond = boundaryConditions[0];
   elements[elemCount - 1].cond = boundaryConditions[1];
}

void NewtonMethod::BuildMatrixNewton()
{
   // Обнуление глобальной матрицы
   for (int i = 0; i < n; i++)
   {
      matrix.di[i] = 0.0;
      for (int j = 0; j < 2; j++)
      {
         matrix.ggl[i][j] = 0.0;
         matrix.ggu[i][j] = 0.0;
      }
   }

   // Сборка из элементов
   for (int i = 0; i < elemCount; i++)
   {
      // Вычисляем линеаризованную локальную матрицу
      elements[i].BuildLocalMatrixNewton(q);

      // Ансамблирование в глобальную матрицу
      for (int j = 0; j < 3; j++)
      {
         int global_j;
         if (i == 0)
            global_j = 3 * i + j;
         else
            global_j = 2 * i + j;

         // Диагональные элементы
         matrix.di[global_j] += elements[i].localMatrix[j][j];

         // Внедиагональные элементы
         if (j == 0 && i > 0)
         {
            // Связь с предыдущим элементом
            matrix.ggl[global_j][1] += elements[i].localMatrix[j][1];
            matrix.ggl[global_j][0] += elements[i].localMatrix[j][2];
         }
         else if (j == 1)
         {
            matrix.ggl[global_j][1] += elements[i].localMatrix[j][0];
            matrix.ggl[global_j][0] += elements[i].localMatrix[j][2];
            matrix.ggu[global_j][1] += elements[i].localMatrix[j][0];
         }
         else if (j == 2 && i < elemCount - 1)
         {
            matrix.ggu[global_j][0] += elements[i].localMatrix[j][0];
            matrix.ggu[global_j][1] += elements[i].localMatrix[j][1];
         }
      }
   }
}

void NewtonMethod::BuildRHSNewton()
{
   // Обнуление глобального вектора правой части
   for (int i = 0; i < n; i++)
   {
      b[i] = 0.0;
   }

   // Сборка из элементов
   for (int i = 0; i < elemCount; i++)
   {
      // Вычисляем локальную правую часть ( -R_i )
      elements[i].BuildLocalBNewton(q);

      // Ансамблирование
      for (int j = 0; j < 3; j++)
      {
         int global_j;
         if (i == 0)
            global_j = 3 * i + j;
         else
            global_j = 2 * i + j;

         b[global_j] += elements[i].localB[j];
      }
   }
}

double NewtonMethod::ComputeResidualNorm()
{
   // Вычисляем невязку R = A(q)*q - b
   // Сначала собираем A(q) и b(q) как в методе простой итерации
   std::vector<double> Aq(n, 0.0);

   // Сборка матрицы A(q) (без линеаризации)
   for (int i = 0; i < elemCount; i++)
   {
      elements[i].BuildLocalMatrix();  // обычная матрица, не линеаризованная

      for (int j = 0; j < 3; j++)
      {
         int global_j;
         if (i == 0)
            global_j = 3 * i + j;
         else
            global_j = 2 * i + j;

         // Умножаем на q
         for (int k = 0; k < 3; k++)
         {
            int global_k;
            if (i == 0)
               global_k = 3 * i + k;
            else
               global_k = 2 * i + k;

            Aq[global_j] += elements[i].localMatrix[j][k] * q[global_k];
         }
      }
   }

   // Сборка вектора b(q)
   std::vector<double> b_res(n, 0.0);
   for (int i = 0; i < elemCount; i++)
   {
      elements[i].BuildLocalB();

      for (int j = 0; j < 3; j++)
      {
         int global_j;
         if (i == 0)
            global_j = 3 * i + j;
         else
            global_j = 2 * i + j;

         b_res[global_j] += elements[i].localB[j];
      }
   }

   // Вычисляем норму невязки
   double norm = 0.0;
   for (int i = 0; i < n; i++)
   {
      double R = Aq[i] - b_res[i];
      norm += R * R;
   }

   return sqrt(norm);
}

void NewtonMethod::SolveNewton()
{
   std::cout << "\n=== Newton Method ===" << std::endl;

   for (int iter = 0; iter < maxIter; iter++)
   {
      newtonIter = iter + 1;

      // 1. Сборка линеаризованной системы
      BuildMatrixNewton();    // A_L
      BuildRHSNewton();       // b_L = -R(q)

      // 2. Применение граничных условий
      ApplyBoundaryConditions();

      // 3. Решение СЛАУ: A_L * dq = b_L
      lu = LU(matrix, b);
      lu.calcLU();
      lu.calcY();
      lu.calcQ();
      dq = lu.q;

      // 4. Обновление решения
      double norm_dq = 0.0;
      for (int i = 0; i < n; i++)
      {
         q[i] += dq[i];
         norm_dq += dq[i] * dq[i];
      }
      norm_dq = sqrt(norm_dq);

      // 5. Вычисление невязки
      double residual = ComputeResidualNorm();

      // 6. Вывод информации
      std::cout << "Iteration " << iter + 1
         << ": |dq| = " << norm_dq
         << ", |R| = " << residual << std::endl;

      // 7. Проверка сходимости
      if (norm_dq < eps && residual < eps)
      {
         std::cout << "\nNewton method converged in " << iter + 1 << " iterations!" << std::endl;
         break;
      }

      if (iter == maxIter - 1)
      {
         std::cout << "\nWarning: Newton method did not converge in " << maxIter << " iterations!" << std::endl;
      }
   }
}