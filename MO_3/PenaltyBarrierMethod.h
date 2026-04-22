#pragma once

#include "Gauss.h"

class PenaltyBarrierMethod
{
public:
   enum class Mode
   {
      Penalty,
      Barrier
   };

private:
   Mode mode;
   int n;
   double* x0;
   double eps;
   double r0;
   double beta;
   int modifierVariant;
   int strategyVariant;
   int maxOuterIter;
   int lastOuterIterCount;
   int totalInnerIterCount;
   long long totalFunctionEvalCount;
   double lastObjectiveValue;
   double lastMeasure;
   double lastCoefficient;
   Gauss gauss;
   Func baseFunc;

   void UpdateCoefficient();

public:
   PenaltyBarrierMethod();
   ~PenaltyBarrierMethod();

   void InitPenalty(int n1, int constraintType);
   void InitBarrier(int n1);
   void SetStartPoint(const double* startPoint);
   void SetTolerance(double epsValue);
   void SetParameters(double startPenalty, double growthFactor);
   void SetModifierVariant(int type);
   void SetStrategyVariant(int type);

   int GetLastOuterIterationCount() const;
   int GetTotalInnerIterationCount() const;
   long long GetTotalFunctionEvaluationCount() const;
   double GetLastObjectiveValue() const;
   double GetLastMeasure() const;
   double GetLastCoefficient() const;

   double* Solver();
};
