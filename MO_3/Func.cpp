#include "Func.h"
#include <algorithm>
#include <cmath>
#include <limits>

void Func::SetFunctionType(int type)
{
   functionType = type;
}

void Func::SetConstraintType(int type)
{
   constraintType = type;
}

void Func::SetModifierType(int type)
{
   modifierType = type;
}

void Func::SetModifierVariant(int type)
{
   modifierVariant = type;
}

void Func::SetCoefficient(double value)
{
   coefficient = value;
}

void Func::ResetEvaluationCount()
{
   evaluationCount = 0;
}

long long Func::GetEvaluationCount() const
{
   return evaluationCount;
}

double Func::Objective(double* x) const
{
   ++evaluationCount;

   switch (functionType)
   {
   case 1:
      return 100.0 * (x[1] - x[0]) * (x[1] - x[0]) + (1.0 - x[0]) * (1.0 - x[0]);
   case 2:
      return 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1.0 - x[0]) * (1.0 - x[0]);
   case 3:
      return -2.0 / (1.0 + (x[0] - 1.0) * (x[0] - 1.0) / 4.0 + (x[1] - 2.0) * (x[1] - 2.0))
         - 1.0 / (1.0 + (x[0] - 3.0) * (x[0] - 3.0) / 9.0 + (x[1] - 1.0) * (x[1] - 1.0) / 9.0);
   case 4:
      return 2.0 * (x[0] - x[1]) * (x[0] - x[1]) + 14.0 * (x[1] - 3.0) * (x[1] - 3.0);
   default:
      return 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1.0 - x[0]) * (1.0 - x[0]);
   }
}

double Func::ConstraintValue(double* x) const
{
   switch (constraintType)
   {
   case 1:
      return x[0] - x[1] + 1.0;
   case 2:
      return x[0] + x[1];
   default:
      return 0.0;
   }
}

double Func::ConstraintViolation(double* x) const
{
   const double value = ConstraintValue(x);

   if (constraintType == 1)
      return std::max(0.0, value);

   if (constraintType == 2)
      return std::fabs(value);

   return 0.0;
}

double Func::func(double* x) const
{
   const double objective = Objective(x);

   if (modifierType == 0 || constraintType == 0)
      return objective;

   if (modifierType == 1)
   {
      const double violation = ConstraintViolation(x);
      switch (modifierVariant)
      {
      case 1:
         return objective + coefficient * violation;
      case 2:
         return objective + coefficient * violation * violation;
      case 3:
         return objective + coefficient * violation * violation * violation * violation;
      default:
         return objective + coefficient * violation * violation;
      }
   }

   if (modifierType == 2)
   {
      const double value = ConstraintValue(x);
      if (constraintType != 1 || value >= 0.0)
         return std::numeric_limits<double>::infinity();

      switch (modifierVariant)
      {
      case 1:
         return objective - coefficient * std::log(-value);
      case 2:
         return objective + coefficient * (-1.0 / value);
      case 3:
         return objective + coefficient * (1.0 / (value * value));
      default:
         return objective - coefficient * std::log(-value);
      }
   }

   return objective;
}

double Func::funcInDirection(double* xk, int coord, double lambda) const
{
   double point[2] = { xk[0], xk[1] };
   point[coord] += lambda;
   return func(point);
}
