#include "Mesh.h"
#include "SimpleIteration.h"
#include "Newton.h"
#include <iostream>
#include "LU.h"

int main()
{
   int method, Id;

   std::cout << "Method (1 - Simple iteration, 2 - Newton): ";
   std::cin >> method;

   std::cout << "Function (1, 2, ...): ";
   std::cin >> Id;

   switch (method)
   {
   case 1:
   {
      SimpleIteration solver;
      solver.functions.id = Id;

      solver.Input();
      solver.Solve();

      for (double v : solver.q)
         std::cout << v << "\t";

      break;
   }

   case 2:
   {
      Newton solver;
      solver.functions.id = Id;

      solver.Input();
      solver.Solve();

      for (double v : solver.q)
         std::cout << v << "\t";

      break;
   }
   }

   return 0;
}