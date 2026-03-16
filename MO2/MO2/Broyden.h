#pragma once
#include "Func.h"
#include "VecOperations.h"
#include "MatrixOperations.h"
class Broyden
{
   int n = 2;
   double *x0;
   double *xk;
   double *xk_1;
   double **etak;
   double lambdak;
   double eps = 1e-7;
   double eps1 = 1e-7;
   Func f;
   VecOperations opsVec;
   MatrixOperations opsMat;
public:
   void Init( int n1 );
   double *Solver( );
   void CalcEtak( );
   double Min( double eps, double an, double bn );
};

