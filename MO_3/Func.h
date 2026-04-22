#pragma once

class Func
{
   int functionType;
   int constraintType;
   int modifierType;
   int modifierVariant;
   double coefficient;
   mutable long long evaluationCount;
public:
   Func()
      : functionType(2), constraintType(0), modifierType(0), modifierVariant(2), coefficient(0.0), evaluationCount(0)
   {
   }

   void SetFunctionType(int type);
   void SetConstraintType(int type);
   void SetModifierType(int type);
   void SetModifierVariant(int type);
   void SetCoefficient(double value);
   void ResetEvaluationCount();
   long long GetEvaluationCount() const;

   double Objective(double* x) const;
   double ConstraintValue(double* x) const;
   double ConstraintViolation(double* x) const;
   double func(double* x) const;
   double funcInDirection(double* xk, int coord, double lambda) const;
};
