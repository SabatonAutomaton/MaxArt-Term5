#include "Mesh.h"
#include "SimpleIteration.h"
#include "Newton.h"
#include <iostream>
#include "LU.h"

int main()
{
   int method;
   std::cout << "Choose method: 1 - Simple Iteration, 2 - Newton\n";
   std::cin >> method;

   switch (method)
   {
   case 1:
   {
      SimpleIteration solver;
      solver.Input();
      solver.Solve();

      for (int i = 0; i < solver.n; i++)
         std::cout << solver.q[i] << "\t";

      break;
   }

   case 2:
   {
      Newton solver;
      solver.Input();
      solver.Solve();

      for (int i = 0; i < solver.n; i++)
         std::cout << solver.q[i] << "\t";

      break;
   }

   default:
      std::cout << "Unknown method\n";
   }

   return 0;
}