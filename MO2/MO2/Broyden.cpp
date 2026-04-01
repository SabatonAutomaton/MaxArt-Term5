#include "Broyden.h"
#include "iostream"
#include <iomanip>
#include <fstream>
#include <cmath>

void Broyden::Init( int n1 )
{
   n = n1;
   x0 = new double[2] { -1, 6 };
   xk = new double[n] { 0 };
   xk_1 = new double[n] { 0 };
   etak = new double *[n];
   for ( int i = 0; i < n; i++ )
   {
      etak[i] = new double[n];
      for ( int j = 0; j < n; j++ )
         etak[i][j] = ( i == j ) ? 1.0 : 0.0;
   }
}

void Broyden::SetFunctionType( int type )
{
   f.SetFunctionType( type );
}

void Broyden::FindInterval( double &a, double &b )
{
   const double H0 = 0.1;
   const double TAU = 2.0;
   const double MIN_STEP = 1e-12;
   const int MAX_ITER = 100;

   a = 0.0;
   double h = H0;
   double f0 = f.funcInDirectionEta( xk, 0.0, etak );
   double f1 = f.funcInDirectionEta( xk, h, etak );

   while ( f1 >= f0 && h > MIN_STEP )
   {
      h *= 0.5;
      f1 = f.funcInDirectionEta( xk, h, etak );
   }

   if ( h <= MIN_STEP )
   {
      b = H0;
      return;
   }

   double lambda_prev = 0.0;
   double lambda_curr = h;
   double f_curr = f1;

   for ( int iter = 0; iter < MAX_ITER; iter++ )
   {
      double lambda_next = lambda_curr * TAU;
      double f_next = f.funcInDirectionEta( xk, lambda_next, etak );

      if ( f_next >= f_curr )
      {
         a = lambda_prev;
         b = lambda_next;
         return;
      }

      lambda_prev = lambda_curr;
      lambda_curr = lambda_next;
      f_curr = f_next;
   }

   b = lambda_curr;
}

double *Broyden::Solver( )
{
   double a = 0, b = 0;
   xk = x0;
   int k = 0;
   std::ofstream file1( "out1.txt" );
   file1 << "i\t" << "x\t" << "y\t" << "f\t" << "lambda\t" << "|xi - xi-1|\t" << "|yi - yi-1|\t" << "|fi - fi-1|\t" << "eta1\teta2\teta3\teta4\t" << "grad1\tgrad2" << std::endl;
   std::ofstream file2( "out2.txt" );
   file2 << 0 << ' ' << xk[0] << ' ' << xk[1] << std::endl;
   while ( opsVec.DotProduct( f.gradFunc( xk ), f.gradFunc( xk ) ) >= eps )
   {
      for ( int i = 0; i < n; i++ ) xk_1[i] = xk[i];
      FindInterval( a, b );
      lambdak = Parabola( eps1, a, b );
      xk = opsVec.AddVec( xk_1, opsVec.MultVecScal( opsMat.MultMatVec( etak, f.gradFunc( xk_1 ) ), -lambdak ) );
      CalcEtak( );
      k++;
      std::cout << k << "\t";
      file1 << k << "\t";
      file2 << k << " ";
      for ( int i = 0; i < n; i++ )
      {
         std::cout << std::setprecision( 16 ) << xk[i] << ' ';
         file1 << std::setprecision( 16 ) << xk[i] << '\t';
         file2 << std::setprecision( 16 ) << xk[i] << ' ';
      }
      std::cout << '\t';
      for ( int i = 0; i < n; i++ )
         for ( int j = 0; j < n; j++ )
            std::cout << std::setprecision( 16 ) << etak[i][j] << ' ';
      file1 << std::setprecision( 16 ) << abs( f.func( xk ) ) << '\t' << lambdak << '\t';
      file1 << abs( xk[0] - xk_1[0] ) << '\t' << abs( xk[1] - xk_1[1] ) << '\t' << abs( f.func( xk ) - f.func( xk_1 ) ) << '\t';
      file1 << etak[0][0] << '\t' << etak[0][1] << '\t' << etak[1][0] << '\t' << etak[1][1] << '\t' << f.gradFunc( xk )[0] << '\t' << f.gradFunc( xk )[1];
      std::cout << std::endl;
      file1 << std::endl;
      file2 << std::endl;
   }
   return xk;
}

void Broyden::CalcEtak( )
{
   double *deltaxk = opsVec.AddVec( xk, opsVec.MultVecScal( xk_1, -1 ) );
   double *deltagk = opsVec.AddVec( f.gradFunc( xk ), opsVec.MultVecScal( f.gradFunc( xk_1 ), -1 ) );
   double *buf = opsVec.AddVec( deltaxk, opsMat.MultMatVec( etak, opsVec.MultVecScal( deltagk, -1 ) ) );
   double denomenator = 1 / opsVec.DotProduct( buf, deltagk );
   etak = opsMat.AddMatrix( opsMat.MultMatScal( opsMat.MultVecVecT( buf, buf ), denomenator ), etak );
   delete[] deltaxk;
   delete[] deltagk;
   delete[] buf;
}

double Broyden::Min( double eps, double an, double bn )
{
   double x, x1, x2;
   double delta = 0.5 * eps;
   int n = 0;
   while ( bn - an >= eps )
   {
      n++;
      x = 0.5 * ( an + bn );
      x1 = x - 0.5 * delta;
      x2 = x + 0.5 * delta;
      double f1 = f.funcInDirectionEta( xk, x1, etak );
      double f2 = f.funcInDirectionEta( xk, x2, etak );
      if ( f1 < f2 )
      {
         bn = x2;
      }
      else if ( f1 > f2 )
      {
         an = x1;
      }
      else
      {
         break;
      }
   }
   return x;
}

double Broyden::GoldenRatio( double eps, double an, double bn )
{
   double x, x1, x2;
   int n = 1;
   double const C = 0.38196601125010515;
   x1 = an + C * ( bn - an );
   x2 = bn - C * ( bn - an );
   double f1 = f.funcInDirectionEta( xk, x1, etak );
   double f2 = f.funcInDirectionEta( xk, x2, etak );
   while ( bn - an >= eps )
   {
      x = 0.5 * ( an + bn );
      n++;
      if ( f1 < f2 )
      {
         bn = x2;
         x2 = x1;
         f2 = f1;
         x1 = an + C * ( bn - an );
         f1 = f.funcInDirectionEta( xk, x1, etak );
      }
      else if ( f1 > f2 )
      {
         an = x1;
         x1 = x2;
         f1 = f2;
         x2 = bn - C * ( bn - an );
         f2 = f.funcInDirectionEta( xk, x2, etak );
      }
      else
      {
         break;
      }
   }
   return x;
}

double Broyden::Parabola( double eps, double an, double bn )
{
   double x1 = an, x2 = ( an + bn ) / 2.0, x3 = bn;
   double x = x2 + 10 * eps;
   double f0, f1 = f.funcInDirectionEta( xk, x1, etak ), f2 = f.funcInDirectionEta( xk, x2, etak ), f3 = f.funcInDirectionEta( xk, x3, etak );
   double numer, denom;
   int k = 0;
   while ( fabs( x3 - x1 ) >= eps && k < maxIter )
   {
      numer = ( x2 - x1 ) * ( x2 - x1 ) * ( f2 - f3 ) - ( x2 - x3 ) * ( x2 - x3 ) * ( f2 - f1 );
      denom = ( x2 - x1 ) * ( f2 - f3 ) - ( x2 - x3 ) * ( f2 - f1 );

      if ( fabs( denom ) < 1e-14 )
         x = x2;
      else
         x = x2 - 0.5 * numer / denom;

      if ( x <= x1 || x >= x3 || !std::isfinite( x ) )
         x = 0.5 * ( x1 + x3 );

      if ( fabs( x - x2 ) < eps )
         break;

      f0 = f.funcInDirectionEta( xk, x, etak );
      if ( x <= x2 )
         if ( f0 <= f2 )
         {
            x3 = x2;
            f3 = f2;
            x2 = x;
            f2 = f0;
         }
         else
         {
            x1 = x;
            f1 = f0;
         }
      else
         if ( f0 <= f2 )
         {
            x1 = x2;
            x2 = x;
            f1 = f2;
            f2 = f0;
         }
         else
         {
            x3 = x;
            f3 = f0;
         }
      k++;
   }
   return x2;
}
