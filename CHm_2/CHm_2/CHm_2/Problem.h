#pragma once
#include "Mesh.h"
#include "SimpleIteration.h"

class Problem
{
   double eps = 1e-16;
   int maxIter = 100;

   //параметры задачи
   double lambda( double u );
   double sigma = 1;
   double *f = nullptr;

   //Краевые условия
   int boundaryConditions[2];
   double theta( double u );
   double beta( double u );
   double uBeta( double u );

   //сетка
   Mesh mesh;

   SimpleIteration simpleIteration;

   void input( );
};

