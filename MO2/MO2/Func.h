#pragma once
#include "VecOperations.h"
#include <iostream>
class Func
{
   VecOperations ops;
public:
   int n=2;
   double func(double *x );
   double *gradFunc(double *x );
   double funcInDirection(double *xk, double lambda, double *sk );
};

