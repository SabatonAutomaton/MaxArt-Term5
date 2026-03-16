#include "Func.h"

double Func::func( double *x )
{
   return 100 * ( x[1] - x[0] ) * ( x[1] - x[0] ) + ( 1 - x[0] ) * ( 1 - x[0] );
   //return 100 * ( x[1] - x[0]* x[0] ) * ( x[1] - x[0]* x[0] ) + ( 1 - x[0] ) * ( 1 - x[0] );
   //return x[0] * x[0] + x[1] * x[1];
}

double *Func::gradFunc( double *x )
{
   double *res = new double[2] {-200 * ( x[1] - x[0] ) - 2 * ( 1 - x[0] ), 200 * ( x[1] - x[0] )};
   //double *res = new double[2] {-400 * x[0]*( x[1] - x[0]*x[0] ) - 2 * ( 1 - x[0] ), 200 * ( x[1] - x[0]*x[0] )};
   //double *res = new double[2] { 2 * x[0], 2 * x[1]};
   return res;
}

double Func::funcInDirection( double *xk, double lambda, double *sk )
{
   return func( ops.AddVec( xk, ops.MultVecScal( sk, lambda ) ) );
}
