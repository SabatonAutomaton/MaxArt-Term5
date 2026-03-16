#pragma once
#include "VecOperations.h"
#include "MatrixOperations.h"
#include <iostream>
class Func
{
   VecOperations opsVec;
   MatrixOperations opsMat;
public:
   int n=2;
   double func(double *x );
   double *gradFunc(double *x );
   double funcInDirection(double *xk, double lambda, double *sk );
   double funcInDirectionEta( double *xk, double lambda, double **etak);
};

