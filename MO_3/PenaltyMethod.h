#pragma once
#include "Gauss.h"

class PenaltyMethod
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
   double lastViolation;
   double lastCoefficient;
   Gauss gauss;
   Func baseFunc;
public:
   PenaltyMethod();
   ~PenaltyMethod();
   void Init(int n1);
   void SetConstraintType(int type);
   void SetStartPoint(const double* startPoint);
   void SetTolerance(double epsValue);
   void SetParameters(double startPenalty, double growthFactor);
   void SetModifierVariant(int type);
   void SetStrategyVariant(int type);
   int GetLastOuterIterationCount() const;
   int GetTotalInnerIterationCount() const;
   long long GetTotalFunctionEvaluationCount() const;
   double GetLastObjectiveValue() const;
   double GetLastViolation() const;
   double GetLastCoefficient() const;
   double* Solver();
};
