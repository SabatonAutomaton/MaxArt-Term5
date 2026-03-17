#pragma once
#include "Func.h"
#include "VecOperations.h"
#include <iostream>
class MSGFR
{
   int n1 = 2;
   int n;
   double *sk;
   double *sk_1;
   double lambdak;
   double *x0;
   double *xk;
   double *xk_1;
   double eps = 1e-7;
   double eps1 = 1e-16;
   Func f;
   VecOperations ops;
public:
   void Init( int n );
   void FindInterval(double& a, double& b);
   double *Solver( );
   double CalcW( );
   double Min( double eps, double an, double bn );
   double GoldenRatio( double eps, double an, double bn );
};

