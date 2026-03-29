#pragma once
#include <cmath>

class Functions
{
public:
   int id = 1;

   double u(double x);
   double f(double x);
   double lambda(double u);
   double dlambda(double u);
};