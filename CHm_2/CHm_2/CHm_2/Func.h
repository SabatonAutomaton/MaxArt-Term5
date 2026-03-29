#pragma once
#include <cmath>

class Functions
{
public:
   int id = 1;
   double currentTime = 0.0;

   double u(double x);
   double f(double x);
   double lambda(double u);
   double dlambda(double u);
   double sigma(double u);
   double sigmaTime(double t);
   double theta(double u);
   double beta(double u);
   double dbeta(double u);
   double uBeta(double u);
   double boundaryValue(int number, double x);
};
