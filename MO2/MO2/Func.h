#pragma once
#include "VecOperations.h"
#include "MatrixOperations.h"
#include <iostream>
class Func
{
   VecOperations opsVec;
   MatrixOperations opsMat;
   int functionType = 2;
public:
   int n=2;
   void SetFunctionType( int type );
   double func(double *x );
   double *gradFunc(double *x );
   double funcInDirection(double *xk, double lambda, double *sk );
   double funcInDirectionEta( double *xk, double lambda, double **etak);
};

