#include "Problem.h"
#include "LU.h"
#include "MSG.h"
#include "SLE.h"
#include "SLEAssistant.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

static double RunLUTest( )
{
   SparseMatrix a;
   a.n = 4;
   a.ig = { 0, 0, 1, 3, 6 };
   a.jg = { 0, 0, 1, 0, 1, 2 };

   a.di = { 10.0, 11.0, 12.0, 13.0 };
   a.ggl = { 4.0, 2.0, 3.0, 1.0, 2.0, 3.0 };
   a.ggu = { 2.0, 3.0, 1.0, 1.0, 2.0, 1.0 };

   const std::vector<double> xTrue = { 1.0, -2.0, 0.5, 3.0 };
   std::vector<std::vector<double>> dense = {
      { 10.0, 2.0, 3.0, 1.0 },
      { 4.0, 11.0, 1.0, 2.0 },
      { 2.0, 3.0, 12.0, 1.0 },
      { 1.0, 2.0, 3.0, 13.0 }
   };

   std::vector<double> b( 4, 0.0 );
   for ( int i = 0; i < 4; ++i )
   {
      for ( int j = 0; j < 4; ++j )
      {
         b[i] += dense[i][j] * xTrue[j];
      }
   }

   LU lu( a );
   std::vector<double> y( 4, 0.0 );
   std::vector<double> x( 4, 0.0 );
   SLEAssistant::LUDirect( lu, b, y );
   SLEAssistant::LUReverse( lu, y, x );

   double maxErr = 0.0;
   for ( int i = 0; i < 4; ++i )
   {
      maxErr = std::max( maxErr, std::abs( x[i] - xTrue[i] ) );
   }

   std::cout << "LU test max|x-x_true| = " << maxErr << "\n";
   return maxErr;
}

static double RunMSGDiagonalTest( )
{
   // Симметричная положительно определённая матрица для теста МСГ.
   SparseMatrix a;
   a.n = 4;
   a.ig = { 0, 0, 1, 3, 6 };
   a.jg = { 0, 0, 1, 0, 1, 2 };

   a.di = { 10.0, 11.0, 12.0, 13.0 };
   a.ggl = { 2.0, 3.0, 1.0, 1.0, 2.0, 1.0 };
   a.ggu = a.ggl;

   const std::vector<double> xTrue = { 1.0, -2.0, 0.5, 3.0 };
   std::vector<std::vector<double>> dense = {
      { 10.0, 2.0, 3.0, 1.0 },
      { 2.0, 11.0, 1.0, 2.0 },
      { 3.0, 1.0, 12.0, 1.0 },
      { 1.0, 2.0, 1.0, 13.0 }
   };

   std::vector<double> b( 4, 0.0 );
   for ( int i = 0; i < 4; ++i )
   {
      for ( int j = 0; j < 4; ++j )
      {
         b[i] += dense[i][j] * xTrue[j];
      }
   }

   SLE sle;
   sle.matrix = a;
   sle.f = b;
   sle.x.assign( 4, 0.0 );
   sle.maxIter = 1000;
   sle.eps = 1e-14;

   const std::vector<double> x = MSG::SolutionWithDiagonalConditioning( sle );

   double maxErr = 0.0;
   for ( int i = 0; i < 4; ++i )
   {
      maxErr = std::max( maxErr, std::abs( x[i] - xTrue[i] ) );
   }

   std::cout << "MSG(diag) test max|x-x_true| = " << maxErr << "\n";
   return maxErr;
}

int main( )
{
   RunLUTest( );
   RunMSGDiagonalTest( );
   Problem problem;
   problem.Run( );
   return 0;
}
