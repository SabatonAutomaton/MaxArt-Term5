#pragma once
#include "Gauss.h"

class BarrierMethod
{
private:
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
   double lastBoundaryGap;
   double lastCoefficient;
   Gauss gauss;
   Func baseFunc;
public:
   BarrierMethod();
   ~BarrierMethod();
   void Init(int n1);
   void SetStartPoint(const double* startPoint);
   void SetTolerance(double epsValue);
   void SetParameters(double startBarrier, double decreaseFactor);
   void SetModifierVariant(int type);
   void SetStrategyVariant(int type);
   int GetLastOuterIterationCount() const;
   int GetTotalInnerIterationCount() const;
   long long GetTotalFunctionEvaluationCount() const;
   double GetLastObjectiveValue() const;
   double GetLastBoundaryGap() const;
   double GetLastCoefficient() const;
   double* Solver();
};
