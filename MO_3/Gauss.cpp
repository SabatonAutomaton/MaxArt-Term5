#include "Gauss.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

Gauss::Gauss()
{
   n = 2;
   xk = nullptr;
   xk_prev = nullptr;
   x0 = nullptr;
   eps_f = 1e-7;
   eps_x = 1e-7;
   maxIter = 10000;
   lastIterCount = 0;
}

Gauss::~Gauss()
{
   delete[] xk;
   delete[] xk_prev;
   delete[] x0;
}

void Gauss::Init(int n1)
{
   n = n1;
   x0 = new double[n];
   xk = new double[n];
   xk_prev = new double[n];

   for (int i = 0; i < n; i++)
   {
      x0[i] = -1.0;
      xk[i] = x0[i];
      xk_prev[i] = x0[i];
   }
}

void Gauss::SetFunctionType(int type)
{
   f.SetFunctionType(type);
}

void Gauss::SetConstraintType(int type)
{
   f.SetConstraintType(type);
}

void Gauss::SetModifierType(int type)
{
   f.SetModifierType(type);
}

void Gauss::SetModifierVariant(int type)
{
   f.SetModifierVariant(type);
}

void Gauss::SetCoefficient(double value)
{
   f.SetCoefficient(value);
}

void Gauss::SetStartPoint(const double* startPoint)
{
   for (int i = 0; i < n; i++)
   {
      x0[i] = startPoint[i];
      xk[i] = startPoint[i];
      xk_prev[i] = startPoint[i];
   }
}

void Gauss::SetTolerances(double epsValue)
{
   eps_f = epsValue;
   eps_x = epsValue;
}

int Gauss::GetLastIterationCount() const
{
   return lastIterCount;
}

long long Gauss::GetFunctionEvaluationCount() const
{
   return f.GetEvaluationCount();
}

double Gauss::GetCurrentFunctionValue() const
{
   return f.func(xk);
}

void Gauss::FindInterval(int coord, double& a, double& b)
{
   const double initialStep = 1.0;
   const double minStep = 1e-12;
   const double expandFactor = 2.0;

   double h = initialStep;
   const double f0 = f.funcInDirection(xk, coord, 0.0);
   double fPlus = f.funcInDirection(xk, coord, h);
   double fMinus = f.funcInDirection(xk, coord, -h);

   while (fPlus >= f0 && fMinus >= f0 && h > minStep)
   {
      h *= 0.5;
      fPlus = f.funcInDirection(xk, coord, h);
      fMinus = f.funcInDirection(xk, coord, -h);
   }

   if (h <= minStep)
   {
      a = -initialStep;
      b = initialStep;
      return;
   }

   const double direction = (fPlus < fMinus) ? 1.0 : -1.0;
   double lambdaPrev = 0.0;
   double lambdaCurr = direction * h;
   double fCurr = f.funcInDirection(xk, coord, lambdaCurr);

   for (int iter = 0; iter < 100; iter++)
   {
      const double lambdaNext = lambdaCurr + direction * h * expandFactor;
      const double fNext = f.funcInDirection(xk, coord, lambdaNext);

      if (fNext >= fCurr)
      {
         a = std::min(lambdaPrev, lambdaNext);
         b = std::max(lambdaPrev, lambdaNext);
         return;
      }

      lambdaPrev = lambdaCurr;
      lambdaCurr = lambdaNext;
      fCurr = fNext;
      h *= expandFactor;
   }

   a = std::min(lambdaPrev, lambdaCurr);
   b = std::max(lambdaPrev, lambdaCurr);
}

double Gauss::GoldenSection(int coord, double left, double right) const
{
   const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
   double x1 = right - (right - left) / phi;
   double x2 = left + (right - left) / phi;
   double f1 = f.funcInDirection(xk, coord, x1);
   double f2 = f.funcInDirection(xk, coord, x2);

   while (std::fabs(right - left) > eps_x)
   {
      if (f1 <= f2)
      {
         right = x2;
         x2 = x1;
         f2 = f1;
         x1 = right - (right - left) / phi;
         f1 = f.funcInDirection(xk, coord, x1);
      }
      else
      {
         left = x1;
         x1 = x2;
         f1 = f2;
         x2 = left + (right - left) / phi;
         f2 = f.funcInDirection(xk, coord, x2);
      }
   }

   return 0.5 * (left + right);
}

double* Gauss::Solver()
{
   std::ofstream file1("gauss_trace.txt");
   std::ofstream file2("out2.txt");
   file1 << "i\tx\ty\tf\tlambda_x\tlambda_y\t|dx|\t|dy|\t|df|" << std::endl;
   file2 << 0 << ' ' << xk[0] << ' ' << xk[1] << std::endl;

   f.ResetEvaluationCount();
   double fPrev = f.func(xk);
   int iter = 0;
   lastIterCount = 0;

   while (iter < maxIter)
   {
      for (int i = 0; i < n; i++)
         xk_prev[i] = xk[i];

      const double fPrevIter = f.func(xk);
      double lambdas[2] = { 0.0, 0.0 };

      for (int coord = 0; coord < n; coord++)
      {
         double a = 0.0;
         double b = 0.0;
         FindInterval(coord, a, b);
         const double lambda = GoldenSection(coord, a, b);
         xk[coord] += lambda;
         lambdas[coord] = lambda;
      }

      const double fCur = f.func(xk);
      ++iter;

      file1 << iter << "\t"
            << std::setprecision(16) << xk[0] << "\t"
            << xk[1] << "\t"
            << fCur << "\t"
            << lambdas[0] << "\t"
            << lambdas[1] << "\t"
            << std::fabs(xk[0] - xk_prev[0]) << "\t"
            << std::fabs(xk[1] - xk_prev[1]) << "\t"
            << std::fabs(fCur - fPrevIter) << std::endl;

      file2 << iter << " "
            << std::setprecision(16) << xk[0] << " "
            << xk[1] << std::endl;

      const bool xSmall = std::fabs(xk[0] - xk_prev[0]) <= eps_x;
      const bool ySmall = std::fabs(xk[1] - xk_prev[1]) <= eps_x;

      if (std::fabs(fCur - fPrev) <= eps_f && xSmall && ySmall)
         break;

      fPrev = fCur;
   }

   lastIterCount = iter;
   return xk;
}
