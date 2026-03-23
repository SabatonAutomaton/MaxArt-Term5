#pragma once
#include "Basis.h"
#include "Mesh.h";
#include "Element.h"
#include "Matrix.h"
#include "LU.h"

class SimpleIteration
{
public:
   double eps = 1e-16;
   int maxIter = 100;
   int elemCount = 0;
   int n = 0;

   //Краевые условия
   int boundaryConditions[2];
   double Condition(int number, double x );
   void Condition1( );

   std::vector<double> q;

   //матрица
   Matrix matrix;

   //вектор правой части
   std::vector<double> b;

   double deltaT;

   std::vector<Element> elements;
   void Input( );
   void BuildMatrix( );
   void BuildB( );
   void SolveIter( );
   double u(double x );
   double f(double x );

   Mesh mesh;
   Basis basis;
   LU lu;
};