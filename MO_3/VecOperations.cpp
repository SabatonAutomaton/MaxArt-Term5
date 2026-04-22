#include "VecOperations.h"

double* VecOperations::MultVecScal(double* vec, double scal)
{
   double* res = new double[n];
   for (int i = 0; i < n; i++)
      res[i] = vec[i] * scal;
   return res;
}

double* VecOperations::AddVec(double* vec1, double* vec2)
{
   double* res = new double[n];
   for (int i = 0; i < n; i++)
      res[i] = vec1[i] + vec2[i];
   return res;
}

double VecOperations::DotProduct(double* vec1, double* vec2)
{
   double res = 0;
   for (int i = 0; i < n; i++)
      res += vec1[i] * vec2[i];
   return res;
}