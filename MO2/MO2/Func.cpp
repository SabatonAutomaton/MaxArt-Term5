#include "Func.h"

double Func::func( double *x )
{
   //return 100 * ( x[1] - x[0] ) * ( x[1] - x[0] ) + ( 1 - x[0] ) * ( 1 - x[0] );
   return 100 * ( x[1] - x[0]* x[0] ) * ( x[1] - x[0]* x[0] ) + ( 1 - x[0] ) * ( 1 - x[0] );
   //return -2 / ( 1 + ( x[0] - 1 ) * ( x[0] - 1 ) / 4.0 + ( x[1] - 1 ) * ( x[1] - 1 ) ) - 1 / ( 1 + ( x[0] - 3 ) * ( x[0] - 3 ) / 9.0 + ( x[1] - 1 ) * ( x[1] - 1 ) / 9.0 );
}

double *Func::gradFunc( double *x )
{
   //double *res = new double[2] {-200 * ( x[1] - x[0] ) - 2 * ( 1 - x[0] ), 200 * ( x[1] - x[0] )};
   double *res = new double[2] {-400 * x[0]*( x[1] - x[0]*x[0] ) - 2 * ( 1 - x[0] ), 200 * ( x[1] - x[0]*x[0] )};
   //double x1 = x[0], x2 = x[1];

   //// Знаменатели (те же, что в исходной функции)
   //double g1 = 1.0 + ( x1 - 1.0 ) * ( x1 - 1.0 ) / 4.0 + ( x2 - 1.0 ) * ( x2 - 1.0 );
   //double g2 = 1.0 + ( ( x1 - 3.0 ) * ( x1 - 3.0 ) + ( x2 - 1.0 ) * ( x2 - 1.0 ) ) / 9.0;

   //// Производные исходной (положительной) функции
   //double dfdx_plus = -( x1 - 1.0 ) / ( g1 * g1 ) - 2.0 * ( x1 - 3.0 ) / ( 9.0 * g2 * g2 );
   //double dfdy_plus = -4.0 * ( x2 - 1.0 ) / ( g1 * g1 ) - 2.0 * ( x2 - 1.0 ) / ( 9.0 * g2 * g2 );

   //// Для отрицательной функции градиент меняет знак
   //double *res = new double[2];
   //res[0] = -dfdx_plus;  // минус от производной исходной
   //res[1] = -dfdy_plus;

   return res;
}

double Func::funcInDirection( double *xk, double lambda, double *sk )
{
   return func( opsVec.AddVec( xk, opsVec.MultVecScal( sk, lambda ) ) );
}

double Func::funcInDirectionEta( double *xk, double lambda, double **etak )
{
   return func( opsVec.AddVec( xk, opsVec.MultVecScal( opsMat.MultMatVec( etak, gradFunc( xk ) ), -lambda ) ) );
}
