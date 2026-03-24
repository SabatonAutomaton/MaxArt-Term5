#pragma once
#include "Basis.h"
#include "Mesh.h"
#include "Element.h"
#include "Matrix.h"
#include "LU.h"

class NewtonMethod
{
public:
   double eps = 1e-8;
   int maxIter = 100;
   int newtonIter = 0;
   int elemCount = 0;
   int n = 0;

   // Краевые условия
   int boundaryConditions[2];
   double Condition(int number, double x);
   void ApplyBoundaryConditions();

   // Решение
   std::vector<double> q;          // текущее приближение
   std::vector<double> q_prev;     // предыдущее приближение (для невязки)
   std::vector<double> dq;         // приращение

   // Матрица и правая часть для линеаризованной системы
   Matrix matrix;
   std::vector<double> b;

   // Элементы
   std::vector<Element> elements;

   // Методы
   void Input();
   void BuildMatrixNewton();
   void BuildRHSNewton();           // правая часть = -R(q)
   void SolveNewton();

   // Вспомогательные функции
   double u(double x);
   double f(double x);
   double ComputeResidualNorm();    // для проверки сходимости

   Mesh mesh;
   Basis basis;
   LU lu;
};