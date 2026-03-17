#include "Broyden.h"
#include "iostream"
#include <iomanip>
#include <fstream>
#include <cmath>
void Broyden::Init( int n1 )
{
   n = n1;
   x0 = new double[2] { 0, -1 };
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

void Broyden::FindInterval( double &a, double &b )
{
   double h = 0.1;
   double lambda0 = 0.5;

   double f0 = f.funcInDirectionEta( xk, lambda0, etak );
   double f1 = f.funcInDirectionEta( xk, lambda0 + h, etak );

   if ( f1 > f0 )
      h = -h;

   double lambda_prev = lambda0;
   double lambda_curr = lambda0 + h;

   int iter = 0;

   while ( iter < 100 )
   {
      double f_prev = f.funcInDirectionEta( xk, lambda_prev, etak );
      double f_curr = f.funcInDirectionEta( xk, lambda_curr, etak );

      if ( f_curr > f_prev )
      {
         a = std::min( lambda_prev, lambda_curr );
         b = std::max( lambda_prev, lambda_curr );
         return;
      }

      lambda_prev = lambda_curr;
      lambda_curr += h;
      h *= 2;
      iter++;
   }

   a = 0;
   b = 1;
}

double *Broyden::Solver( )
{
   double a = 0, b = 0;
   xk = x0;
   int k = 0;
   const double PI = 3.14159265359;
   std::ofstream file1( "out1.txt" );
   file1 << "i\t" << "x\t" << "y\t" << "f\t" << "lambda\t" << "|xi - xi-1|\t" << "|yi - yi-1|\t" << "|fi - fi-1|\t" << "eta1\teta2\teta3\teta4\t" << "grad1\tgrad2" << std::endl;
   std::ofstream file2( "out2.txt" );
   file2 << 0 << ' ' << xk[0] << ' ' << xk[1] << std::endl;
   while ( opsVec.DotProduct( f.gradFunc( xk ), f.gradFunc( xk ) ) >= eps )
   {
      for ( int i = 0; i < n; i++ ) xk_1[i] = xk[i];
      FindInterval( a, b );
      lambdak = Min( eps1, a, b );
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