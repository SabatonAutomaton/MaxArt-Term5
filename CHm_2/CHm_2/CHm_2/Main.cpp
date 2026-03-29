#include "Mesh.h"
#include "SimpleIteration.h"
#include "Newton.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include "LU.h"

std::string ToCommaString( double value, int precision, bool scientific = false )
{
   std::ostringstream stream;
   if ( scientific )
      stream << std::uppercase << std::scientific;
   else
      stream << std::fixed;

   stream << std::setprecision( precision ) << value;
   std::string text = stream.str( );
   for ( char& ch : text )
   {
      if ( ch == '.' )
         ch = ',';
   }
   return text;
}

template <typename Solver>
void WriteSolutionToFile( const Solver& solver, const std::string& fileName )
{
   std::ofstream file( fileName );
   if ( !file )
      return;

   file << "x\tq\tq*\tq - q*" << std::endl;

   Functions exactFunctions = solver.functions;

   for ( int i = 0; i < solver.n; i++ )
   {
      double x = 0.0;
      if ( i % 2 == 0 )
      {
         x = solver.mesh.meshX[i / 2];
      }
      else
      {
         int left = i / 2;
         x = ( solver.mesh.meshX[left] + solver.mesh.meshX[left + 1] ) / 2.0;
      }

      double exact = exactFunctions.u( x );
      double diff = solver.q[i] - exact;

      file
         << ToCommaString( x, 16 ) << "\t"
         << ToCommaString( solver.q[i], 16 ) << "\t"
         << ToCommaString( exact, 16 ) << "\t"
         << ToCommaString( diff, 2, true ) << std::endl;
   }
}

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
      solver.useTimeDependent = true;

      solver.Input();
      solver.Solve();
      WriteSolutionToFile( solver, "out.txt" );

      for (double v : solver.q)
         std::cout << v << "\t";

      break;
   }

   case 2:
   {
      Newton solver;
      solver.functions.id = Id;
      solver.useTimeDependent = true;

      solver.Input();
      solver.Solve();
      WriteSolutionToFile( solver, "out.txt" );

      for (double v : solver.q)
         std::cout << v << "\t";

      break;
   }
   }

   return 0;
}
