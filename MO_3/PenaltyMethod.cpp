#include "PenaltyMethod.h"
#include <fstream>
#include <iomanip>

PenaltyMethod::PenaltyMethod()
{
   n = 2;
   x0 = nullptr;
   eps = 1e-5;
   r0 = 1.0;
   beta = 10.0;
   modifierVariant = 2;
   strategyVariant = 1;
   maxOuterIter = 50;
   lastOuterIterCount = 0;
   totalInnerIterCount = 0;
   totalFunctionEvalCount = 0;
   lastObjectiveValue = 0.0;
   lastViolation = 0.0;
   lastCoefficient = 0.0;
}

PenaltyMethod::~PenaltyMethod()
{
   delete[] x0;
}

void PenaltyMethod::Init(int n1)
{
   n = n1;
   x0 = new double[n];

   for (int i = 0; i < n; i++)
      x0[i] = -1.0;

   gauss.Init(n);
   gauss.SetFunctionType(4);
   gauss.SetModifierType(1);
   gauss.SetModifierVariant(modifierVariant);
   baseFunc.SetFunctionType(4);
}

void PenaltyMethod::SetConstraintType(int type)
{
   gauss.SetConstraintType(type);
   baseFunc.SetConstraintType(type);
}

void PenaltyMethod::SetStartPoint(const double* startPoint)
{
   for (int i = 0; i < n; i++)
      x0[i] = startPoint[i];
}

void PenaltyMethod::SetTolerance(double epsValue)
{
   eps = epsValue;
   gauss.SetTolerances(epsValue);
}

void PenaltyMethod::SetParameters(double startPenalty, double growthFactor)
{
   r0 = startPenalty;
   beta = growthFactor;
}

void PenaltyMethod::SetModifierVariant(int type)
{
   modifierVariant = type;
   gauss.SetModifierVariant(type);
}

void PenaltyMethod::SetStrategyVariant(int type)
{
   strategyVariant = type;
}

int PenaltyMethod::GetLastOuterIterationCount() const
{
   return lastOuterIterCount;
}

int PenaltyMethod::GetTotalInnerIterationCount() const
{
   return totalInnerIterCount;
}

long long PenaltyMethod::GetTotalFunctionEvaluationCount() const
{
   return totalFunctionEvalCount;
}

double PenaltyMethod::GetLastObjectiveValue() const
{
   return lastObjectiveValue;
}

double PenaltyMethod::GetLastViolation() const
{
   return lastViolation;
}

double PenaltyMethod::GetLastCoefficient() const
{
   return lastCoefficient;
}

double* PenaltyMethod::Solver()
{
   std::ofstream file("out_penalty.txt");
   std::ofstream traj("out_penalty_traj.txt");
   file << "outer\tr\tx\ty\tf\tviolation\tinner_iters" << std::endl;
   traj << 0 << ' ' << x0[0] << ' ' << x0[1] << std::endl;

   gauss.SetStartPoint(x0);
   totalInnerIterCount = 0;
   totalFunctionEvalCount = 0;
   lastOuterIterCount = 0;
   lastCoefficient = r0;

   for (int outerIter = 1; outerIter <= maxOuterIter; outerIter++)
   {
      gauss.SetCoefficient(lastCoefficient);
      double* current = gauss.Solver();
      totalFunctionEvalCount += gauss.GetFunctionEvaluationCount();

      for (int i = 0; i < n; i++)
         x0[i] = current[i];

      lastObjectiveValue = baseFunc.Objective(current);
      lastViolation = baseFunc.ConstraintViolation(current);
      totalInnerIterCount += gauss.GetLastIterationCount();
      lastOuterIterCount = outerIter;

      file << outerIter << "\t"
           << std::setprecision(16) << lastCoefficient << "\t"
           << current[0] << "\t"
           << current[1] << "\t"
           << lastObjectiveValue << "\t"
           << lastViolation << "\t"
           << gauss.GetLastIterationCount() << std::endl;
      traj << outerIter << ' ' << std::setprecision(16) << current[0] << ' ' << current[1] << std::endl;

      if (lastViolation <= eps)
         return current;

      gauss.SetStartPoint(current);
      switch (strategyVariant)
      {
      case 1:
         lastCoefficient *= beta;
         break;
      case 2:
         lastCoefficient += beta;
         break;
      case 3:
         lastCoefficient += 1.0;
         break;
      default:
         lastCoefficient *= beta;
         break;
      }
   }

   double* current = gauss.Solver();
   totalFunctionEvalCount += gauss.GetFunctionEvaluationCount();
   return current;
}
