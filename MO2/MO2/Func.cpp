#include "Func.h"

void Func::SetFunctionType( int type )
{
   functionType = type;
}

double Func::func(double* x)
{
   switch ( functionType )
   {
   case 1:
      return 100.0 * ( x[1] - x[0] ) * ( x[1] - x[0] ) + ( 1.0 - x[0] ) * ( 1.0 - x[0] );
   case 2:
      return 100.0 * ( x[1] - x[0] * x[0] ) * ( x[1] - x[0] * x[0] ) + ( 1.0 - x[0] ) * ( 1.0 - x[0] );
   case 3:
      return - 2.0 / ( 1.0 + ( x[0] - 1.0 ) * ( x[0] - 1.0 ) / 4.0 + ( x[1] - 2.0 ) * ( x[1] - 2.0 ) ) -
         1.0 / ( 1.0 + ( x[0] - 3.0 ) * ( x[0] - 3.0 ) / 9.0 + ( x[1] - 1.0 ) * ( x[1] - 1.0 ) / 9.0 );
   default:
      return 100.0 * ( x[1] - x[0] * x[0] ) * ( x[1] - x[0] * x[0] ) + ( 1.0 - x[0] ) * ( 1.0 - x[0] );
   }
}

double* Func::gradFunc(double* x)
{
   double* res = new double[2];

   switch ( functionType )
   {
   case 1:
      res[0] = -200.0 * ( x[1] - x[0] ) - 2.0 * ( 1.0 - x[0] );
      res[1] = 200.0 * ( x[1] - x[0] );
      break;
   case 2:
      res[0] = -400.0 * x[0] * ( x[1] - x[0] * x[0] ) - 2.0 * ( 1.0 - x[0] );
      res[1] = 200.0 * ( x[1] - x[0] * x[0] );
      break;
   case 3:
   {
      double g1 = 1.0 + ( x[0] - 1.0 ) * ( x[0] - 1.0 ) / 4.0 + ( x[1] - 2.0 ) * ( x[1] - 2.0 );
      double g2 = 1.0 + ( x[0] - 3.0 ) * ( x[0] - 3.0 ) / 9.0 + ( x[1] - 1.0 ) * ( x[1] - 1.0 ) / 9.0;

      res[0] = ( x[0] - 1.0 ) / ( g1 * g1 ) + 2.0 * ( x[0] - 3.0 ) / ( 9.0 * g2 * g2 );
      res[1] = 4.0 * ( x[1] - 2.0 ) / ( g1 * g1 ) + 2.0 * ( x[1] - 1.0 ) / ( 9.0 * g2 * g2 );
      break;
   }
   default:
      res[0] = -400.0 * x[0] * ( x[1] - x[0] * x[0] ) - 2.0 * ( 1.0 - x[0] );
      res[1] = 200.0 * ( x[1] - x[0] * x[0] );
      break;
   }

   return res;
}

double Func::funcInDirection(double* xk, double lambda, double* sk)
{
   return func( opsVec.AddVec( xk, opsVec.MultVecScal( sk, lambda ) ) );
}

double Func::funcInDirectionEta( double *xk, double lambda, double **etak )
{
   return func( opsVec.AddVec( xk, opsVec.MultVecScal( opsMat.MultMatVec( etak, gradFunc( xk ) ), -lambda ) ) );
}
