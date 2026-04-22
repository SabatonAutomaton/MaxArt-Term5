#include "PenaltyBarrierMethod.h"

#include <fstream>
#include <iomanip>
#include <limits>

PenaltyBarrierMethod::PenaltyBarrierMethod()
{
   mode = Mode::Penalty;
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
   lastMeasure = 0.0;
   lastCoefficient = 0.0;
}

PenaltyBarrierMethod::~PenaltyBarrierMethod()
{
   delete[] x0;
}

void PenaltyBarrierMethod::InitPenalty(int n1, int constraintType)
{
   mode = Mode::Penalty;
   n = n1;
   delete[] x0;
   x0 = new double[n];
   for (int i = 0; i < n; ++i)
      x0[i] = -1.0;

   gauss.Init(n);
   gauss.SetFunctionType(4);
   gauss.SetConstraintType(constraintType);
   gauss.SetModifierType(1);
   gauss.SetModifierVariant(modifierVariant);

   baseFunc.SetFunctionType(4);
   baseFunc.SetConstraintType(constraintType);
}

void PenaltyBarrierMethod::InitBarrier(int n1)
{
   mode = Mode::Barrier;
   n = n1;
   delete[] x0;
   x0 = new double[n];
   for (int i = 0; i < n; ++i)
      x0[i] = -1.0;

   gauss.Init(n);
   gauss.SetFunctionType(4);
   gauss.SetConstraintType(1);
   gauss.SetModifierType(2);
   gauss.SetModifierVariant(modifierVariant);

   baseFunc.SetFunctionType(4);
   baseFunc.SetConstraintType(1);
}

void PenaltyBarrierMethod::SetStartPoint(const double* startPoint)
{
   for (int i = 0; i < n; ++i)
      x0[i] = startPoint[i];
}

void PenaltyBarrierMethod::SetTolerance(double epsValue)
{
   eps = epsValue;
   gauss.SetTolerances(epsValue);
}

void PenaltyBarrierMethod::SetParameters(double startPenalty, double growthFactor)
{
   r0 = startPenalty;
   beta = growthFactor;
}

void PenaltyBarrierMethod::SetModifierVariant(int type)
{
   modifierVariant = type;
   gauss.SetModifierVariant(type);
}

void PenaltyBarrierMethod::SetStrategyVariant(int type)
{
   strategyVariant = type;
}

int PenaltyBarrierMethod::GetLastOuterIterationCount() const
{
   return lastOuterIterCount;
}

int PenaltyBarrierMethod::GetTotalInnerIterationCount() const
{
   return totalInnerIterCount;
}

long long PenaltyBarrierMethod::GetTotalFunctionEvaluationCount() const
{
   return totalFunctionEvalCount;
}

double PenaltyBarrierMethod::GetLastObjectiveValue() const
{
   return lastObjectiveValue;
}

double PenaltyBarrierMethod::GetLastMeasure() const
{
   return lastMeasure;
}

double PenaltyBarrierMethod::GetLastCoefficient() const
{
   return lastCoefficient;
}

void PenaltyBarrierMethod::UpdateCoefficient()
{
   if (mode == Mode::Penalty)
   {
      switch (strategyVariant)
      {
      case 1: lastCoefficient *= beta; break;
      case 2: lastCoefficient += beta; break;
      case 3: lastCoefficient += 1.0; break;
      default: lastCoefficient *= beta; break;
      }
      return;
   }

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

double* PenaltyBarrierMethod::Solver()
{
   const bool penaltyMode = (mode == Mode::Penalty);
   std::ofstream file(penaltyMode ? "out_penalty.txt" : "out_barrier.txt");
   std::ofstream traj(penaltyMode ? "out_penalty_traj.txt" : "out_barrier_traj.txt");

   if (penaltyMode)
      file << "outer\tr\tx\ty\tf\tviolation\tinner_iters" << std::endl;
   else
      file << "outer\tr\tx\ty\tf\tboundary_gap\tinner_iters" << std::endl;

   traj << 0 << ' ' << x0[0] << ' ' << x0[1] << std::endl;

   gauss.SetStartPoint(x0);
   baseFunc.ResetEvaluationCount();
   totalInnerIterCount = 0;
   totalFunctionEvalCount = 0;
   lastOuterIterCount = 0;
   lastCoefficient = r0;

   for (int outerIter = 1; outerIter <= maxOuterIter; ++outerIter)
   {
      gauss.SetCoefficient(lastCoefficient);
      double* current = gauss.Solver();
      totalFunctionEvalCount += gauss.GetFunctionEvaluationCount();

      for (int i = 0; i < n; ++i)
         x0[i] = current[i];

      lastObjectiveValue = baseFunc.Objective(current);
      ++totalFunctionEvalCount;
      lastMeasure = penaltyMode ? baseFunc.ConstraintViolation(current) : -baseFunc.ConstraintValue(current);
      totalInnerIterCount += gauss.GetLastIterationCount();
      lastOuterIterCount = outerIter;

      file << outerIter << "\t"
         << std::setprecision(16) << lastCoefficient << "\t"
         << current[0] << "\t"
         << current[1] << "\t"
         << lastObjectiveValue << "\t"
         << lastMeasure << "\t"
         << gauss.GetLastIterationCount() << std::endl;
      traj << outerIter << ' ' << std::setprecision(16) << current[0] << ' ' << current[1] << std::endl;

      if (lastMeasure <= eps)
         return current;

      gauss.SetStartPoint(current);
      UpdateCoefficient();
   }

   double* current = gauss.Solver();
   totalFunctionEvalCount += gauss.GetFunctionEvaluationCount();
   return current;
}
