#include "Func.h"

double Functions::u(double x)
{
   switch (id)
   {
   case 1:
      return x;

   case 2:
      return x * x;

   default:
      return x;
   }
}

double Functions::f(double x)
{
   switch (id)
   {
   case 1:
      return x - 1;

   case 2:
      return 2 * x;

   default:
      return x - 1;
   }
}

double Functions::lambda(double u)
{
   switch (id)
   {
   case 1:
      return 1 + u;

   case 2:
      return 1 + u * u;

   default:
      return 1 + u;
   }
}

double Functions::dlambda(double u)
{
   switch (id)
   {
   case 1:
      return 1;

   case 2:
      return 2 * u;

   default:
      return 1;
   }
}