#pragma once
#include "Func.h"

class Gauss
{
private:
   int n;
   double* xk;
   double* xk_prev;
   double* x0;
   double eps_f;
   double eps_x;
   int maxIter;
   int lastIterCount;
   Func f;

   double GoldenSection(int coord, double left, double right) const;
   void FindInterval(int coord, double& a, double& b);
public:
   Gauss();
   ~Gauss();
   void Init(int n1);
   void SetFunctionType(int type);
   void SetConstraintType(int type);
   void SetModifierType(int type);
   void SetModifierVariant(int type);
   void SetCoefficient(double value);
   void SetStartPoint(const double* startPoint);
   void SetTolerances(double epsValue);
   int GetLastIterationCount() const;
   long long GetFunctionEvaluationCount() const;
   double GetCurrentFunctionValue() const;
   double* Solver();
};
