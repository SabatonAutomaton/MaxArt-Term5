#include "BarrierMethod.h"
#include <fstream>
#include <iomanip>

BarrierMethod::BarrierMethod()
{
   n = 2;
   x0 = nullptr;
   eps = 1e-5;
   r0 = 1.0;
   beta = 10.0;
   modifierVariant = 1;
   strategyVariant = 1;
   maxOuterIter = 50;
   lastOuterIterCount = 0;
   totalInnerIterCount = 0;
   totalFunctionEvalCount = 0;
   lastObjectiveValue = 0.0;
   lastBoundaryGap = 0.0;
   lastCoefficient = 0.0;
}

BarrierMethod::~BarrierMethod()
{
   delete[] x0;
}

void BarrierMethod::Init(int n1)
{
   n = n1;
   x0 = new double[n];

   for (int i = 0; i < n; i++)
      x0[i] = -1.0;

   gauss.Init(n);
   gauss.SetFunctionType(4);
   gauss.SetConstraintType(1);
   gauss.SetModifierType(2);
   gauss.SetModifierVariant(modifierVariant);

   baseFunc.SetFunctionType(4);
   baseFunc.SetConstraintType(1);
}

void BarrierMethod::SetStartPoint(const double* startPoint)
{
   for (int i = 0; i < n; i++)
      x0[i] = startPoint[i];
}

void BarrierMethod::SetTolerance(double epsValue)
{
   eps = epsValue;
   gauss.SetTolerances(epsValue);
}

void BarrierMethod::SetParameters(double startBarrier, double decreaseFactor)
{
   r0 = startBarrier;
   beta = decreaseFactor;
}

void BarrierMethod::SetModifierVariant(int type)
{
   modifierVariant = type;
   gauss.SetModifierVariant(type);
}

void BarrierMethod::SetStrategyVariant(int type)
{
   strategyVariant = type;
}

int BarrierMethod::GetLastOuterIterationCount() const
{
   return lastOuterIterCount;
}

int BarrierMethod::GetTotalInnerIterationCount() const
{
   return totalInnerIterCount;
}

long long BarrierMethod::GetTotalFunctionEvaluationCount() const
{
   return totalFunctionEvalCount;
}

double BarrierMethod::GetLastObjectiveValue() const
{
   return lastObjectiveValue;
}

double BarrierMethod::GetLastBoundaryGap() const
{
   return lastBoundaryGap;
}

double BarrierMethod::GetLastCoefficient() const
{
   return lastCoefficient;
}

double* BarrierMethod::Solver()
{
   std::ofstream file("out_barrier.txt");
   std::ofstream traj("out_barrier_traj.txt");
   file << "outer\tr\tx\ty\tf\tboundary_gap\tinner_iters" << std::endl;
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
      lastBoundaryGap = -baseFunc.ConstraintValue(current);
      totalInnerIterCount += gauss.GetLastIterationCount();
      lastOuterIterCount = outerIter;

      file << outerIter << "\t"
           << std::setprecision(16) << lastCoefficient << "\t"
           << current[0] << "\t"
           << current[1] << "\t"
           << lastObjectiveValue << "\t"
           << lastBoundaryGap << "\t"
           << gauss.GetLastIterationCount() << std::endl;
      traj << outerIter << ' ' << std::setprecision(16) << current[0] << ' ' << current[1] << std::endl;

      if (lastBoundaryGap <= eps)
         return current;

      gauss.SetStartPoint(current);
      switch (strategyVariant)
      {
      case 1:
         lastCoefficient /= beta;
         break;
      case 2:
         lastCoefficient -= beta;
         if (lastCoefficient <= 1e-12)
            lastCoefficient = 1e-12;
         break;
      case 3:
         lastCoefficient /= 2.0;
         break;
      default:
         lastCoefficient /= beta;
         break;
      }
   }

   double* current = gauss.Solver();
   totalFunctionEvalCount += gauss.GetFunctionEvaluationCount();
   return current;
}
