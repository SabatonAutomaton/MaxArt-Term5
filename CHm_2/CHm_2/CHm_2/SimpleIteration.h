#pragma once
#include "Basis.h"
#include "Mesh.h";
#include "Element.h"
#include "Matrix.h"
#include "LU.h"
#include "MatrixOperations.h"
#include "Func.h"

class SimpleIteration
{
public:
   double eps = 1e-16;
   int maxIter = 100;
   int elemCount = 0;
   int n = 0;

   Functions functions;

   //Краевые условия
   int boundaryConditions[2];
   double Condition(int number, double x );
   void Condition1( );

   std::vector<double> q;
   std::vector<double> qOLd;

   //матрица
   Matrix matrix;

   //вектор правой части
   std::vector<double> b;
   std::vector<double> b0;

   double deltaT;

   double w = 0.1;

   std::vector<Element> elements;
   void Input( );
   void BuildMatrix( );
   void BuildB( );
   void SolveIter( );
   void Solve( );
   double CalcResidual( );
   std::vector<double> lambdaNodes;

   Mesh mesh;
   Basis basis;
   LU lu;
   MatrixOperations ops;
};
