#include "MSGFR.h"
#include "iostream"
#include <iomanip>
#include <fstream>
#include <cmath>

void MSGFR::Init( int n1 )
{
   n = n1;
   x0 = new double[2] { 5, 6 };
   sk = ops.MultVecScal( f.gradFunc( x0 ), -1 );
   sk_1 = new double[n];
   xk = new double[n];
   xk_1 = new double[n];
   xk_old = new double[n];
   for ( int i = 0; i < n; i++ )
   {
      xk[i] = x0[i];
      xk_1[i] = x0[i];
      xk_old[i] = x0[i];
   }
}

//void MSGFR::FindInterval( double &a, double &b, double *x )
//{
//   double h = 0.1;
//   double lambda = 0.5;
//
//   double f_current = f.funcInDirection( x, lambda, sk );
//   double f_next = f.funcInDirection( x, lambda + h, sk );
//
//   // Определяем направление убывания
//   if ( f_next > f_current ) {
//      h = -h;
//      f_next = f.funcInDirection( x, lambda + h, sk );
//   }
//
//   // Поиск интервала, где функция начинает возрастать
//   while ( f_next < f_current ) {
//      lambda += h;
//      f_current = f_next;
//      h *= 2.0;
//      f_next = f.funcInDirection( x, lambda + h, sk );
//
//      // Защита от зацикливания
//      if ( fabs( h ) > 100 ) break;
//   }
//
//   // Устанавливаем границы интервала
//   if ( h > 0 ) {
//      a = lambda;
//      b = lambda + h;
//   }
//   else {
//      a = lambda + h;
//      b = lambda;
//   }
//}

void MSGFR::FindInterval( double &a, double &b )
{
   const double H0 = 0.1;        // Начальный шаг
   const double TAU = 2.0;        // Множитель увеличения шага
   const int MAX_ITER = 100;      // Максимальное число итераций

   double lambda = 0.0;           // Текущая точка
   double h = H0;                 // Текущий шаг

   double f_prev = f.funcInDirection( xk, lambda - h, sk );
   double f_curr = f.funcInDirection( xk, lambda, sk );
   double f_next = f.funcInDirection( xk, lambda + h, sk );

   if ( f_prev < f_curr ) 
   {
      h = -h;
   }
   else if ( f_next < f_curr ) { }
   else 
   {
      a = lambda - h;
      b = lambda + h;
      return;
   }

   int iter = 0;
   double lambda_prev = lambda;
   double lambda_curr = lambda + h;
   double f_curr_val, f_next_val;

   while ( iter < MAX_ITER ) {
      f_curr_val = f.funcInDirection( xk, lambda_curr, sk );
      f_next_val = f.funcInDirection( xk, lambda_curr + h, sk );


      if ( f_curr_val < f_next_val ) {
         a = std::min( lambda_prev, lambda_curr + h );
         b = std::max( lambda_prev, lambda_curr + h );

         return;
      }
      lambda_prev = lambda_curr;
      lambda_curr = lambda_curr + h;
      h *= TAU;

      iter++;
   }
   a = 0.0;
   b = 10.0;
}

double *MSGFR::Solver( )
{
   const double PI = 3.14159265359;
   std::ofstream file1( "out1.txt" );
   file1 << "i\t" << "x\t" << "y\t" << "f\t" << "s1\t" << "s2\t" << "lambda\t" << "|xi - xi-1|\t" << "|yi - yi-1|\t" << "|fi - fi-1|\t" << "angle\t" << "grad1\tgrad2" << std::endl;
   std::ofstream file2( "out2.txt" );
   file2 << 0 << ' ' << x0[0] << ' ' << x0[1] << std::endl;
   double a = -1, b = 1;
   //FindInterval( a, b, x0 );
   //lambdak = Parabola( eps1, a, b );
   //xk = ops.AddVec( x0, ops.MultVecScal( sk, lambdak ) );
   int k = 0;
   while ( ops.DotProduct( sk, sk ) >= eps )
   {
      if ( (k+1) % n1 == 0 )
      {
         sk = ops.MultVecScal( f.gradFunc( xk ), -1 );
      }
      FindInterval( a, b);
      lambdak = Parabola( eps1, a, b );
      for ( int i = 0; i < n; i++ ) xk_old[i] = xk[i];
      xk = ops.AddVec( xk, ops.MultVecScal( sk, lambdak ) );
      for ( int i = 0; i < n; i++ ) xk_1[i] = xk_old[i];
      for ( int i = 0; i < n; i++ ) sk_1[i] = sk[i];
      double wk = CalcW( );
      sk = ops.AddVec( ops.MultVecScal( f.gradFunc( xk ), -1 ), ops.MultVecScal( sk, wk ) );
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
      file1 << std::setprecision( 16 ) << f.func( xk ) << '\t' << sk[0] << '\t' << sk[1] << '\t' << lambdak << '\t';
      file1 << abs( xk[0] - xk_1[0] ) << '\t' << abs( xk[1] - xk_1[1] ) << '\t' << abs( f.func( xk ) - f.func( xk_1 ) ) << '\t';
      file1 << std::acos( ( xk[0] * sk[0] + xk[1] * sk[1] ) / ( sqrt( ops.DotProduct( xk, xk ) ) * sqrt( ops.DotProduct( sk, sk ) ) ) ) * 180 / PI << '\t';
      file1 << f.gradFunc( xk )[0] << '\t' << f.gradFunc( xk )[1];
      std::cout << std::endl;
      file1 << std::endl;
      file2 << std::endl;
   }
   return xk;
}

double MSGFR::CalcW( )
{
   double a = ops.DotProduct( f.gradFunc( xk ), f.gradFunc( xk ) );
   double b = ops.DotProduct( f.gradFunc( xk_1 ), f.gradFunc( xk_1 ) );
   return a / b;
}

double MSGFR::Min( double eps, double an, double bn )
{
   double x, x1, x2;
   double delta = 0.5 * eps1;
   int n = 0;
   while ( bn - an >= eps )
   {
      n++;
      x = 0.5 * ( an + bn );
      x1 = x - 0.5 * delta;
      x2 = x + 0.5 * delta;
      double f1 = f.funcInDirection( xk, x1, sk );
      double f2 = f.funcInDirection( xk, x2, sk );
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


double MSGFR::GoldenRatio( double eps, double an, double bn )
{
   double x, x1, x2;
   int n = 1;
   double const C = 0.38196601125010515;
   x1 = an + C * ( bn - an );
   x2 = bn - C * ( bn - an );
   double f1 = f.funcInDirection( xk, x1, sk );
   double f2 = f.funcInDirection( xk, x2, sk );
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
         f1 = f.funcInDirection( xk, x1, sk );
      }
      else if ( f1 > f2 )
      {
         an = x1;
         x1 = x2;
         f1 = f2;
         x2 = bn - C * ( bn - an );
         f2 = f.funcInDirection( xk, x2, sk );
      }
      else
      {
         break;
      }
   }
   return x;
}

double MSGFR::Parabola( double eps, double an, double bn )
{
   double x1 = an, x2 = ( an + bn ) / 2.0, x3 = bn;
   double x = x2 + 10 * eps;
   double f0, f1 = f.funcInDirection( xk, x1, sk ), f2 = f.funcInDirection( xk, x2, sk ), f3 = f.funcInDirection( xk, x3, sk );
   double numer, denom; // числитель, знаменатель
   int k = 0;
   while ( fabs( x2 - x ) >= eps )
   {
      numer = ( x2 - x1 ) * ( x2 - x1 ) * ( f2 - f3 ) - ( x2 - x3 ) * ( x2 - x3 ) * ( f2 - f1 );
      denom = ( x2 - x1 ) * ( f2 - f3 ) - ( x2 - x3 ) * ( f2 - f1 );
      if ( denom == 0 )
         x = x2;
      else
         x = x2 - 0.5 * numer / denom;
      f0 = f.funcInDirection( xk, x, sk );
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
   std::cout << k << std::endl;
   return x2;
}


