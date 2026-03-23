#include "Element.h"

double Element::lambda( double u )
{
   return 1;
}

double Element::theta( double u )
{
   return 0.0;
}

double Element::beta( double u )
{
   return 0.0;
}

double Element::uBeta( double u )
{
   return 0.0;
}

double Element::uh( double x )
{
   return q[0] * basis.basis1( x ) + q[1] * basis.basis2( x ) + q[2] * basis.basis3( x );
}

void Element::BuildLocalMatrix( )
{
   double **quadrature = basis.Quadrature( );
   double jacobian = coords[1] - coords[0];
   for ( int i = 0; i < 3; i++ )
   {
      double wj = quadrature[i][1] * jacobian;
      double psi1 = basis.basis1( quadrature[i][0] );
      double psi2 = basis.basis2( quadrature[i][0] );
      double psi3 = basis.basis3( quadrature[i][0] );
      double gradPsi1 = basis.basis1Grad( quadrature[i][0] ) / jacobian;
      double gradPsi2 = basis.basis2Grad( quadrature[i][0] ) / jacobian;
      double gradPsi3 = basis.basis3Grad( quadrature[i][0] ) / jacobian;
      double lambdas = lambda( uh( quadrature[i][0] ) );
      localMatrix[0][0] += wj * ( lambdas * gradPsi1 * gradPsi1 + sigma * psi1 * psi1 + beta( uh( quadrature[i][0] ) ) * psi1 * psi1 );
      localMatrix[0][1] += wj * ( lambdas * gradPsi1 * gradPsi2 + sigma * psi1 * psi2 + beta( uh( quadrature[i][0] ) ) * psi1 * psi2 );
      localMatrix[0][2] += wj * ( lambdas * gradPsi1 * gradPsi3 + sigma * psi1 * psi3 + beta( uh( quadrature[i][0] ) ) * psi1 * psi3 );
      localMatrix[1][0] += wj * ( lambdas * gradPsi2 * gradPsi1 + sigma * psi2 * psi1 + beta( uh( quadrature[i][0] ) ) * psi2 * psi1 );
      localMatrix[1][1] += wj * ( lambdas * gradPsi2 * gradPsi2 + sigma * psi2 * psi2 + beta( uh( quadrature[i][0] ) ) * psi2 * psi2 );
      localMatrix[1][2] += wj * ( lambdas * gradPsi2 * gradPsi3 + sigma * psi2 * psi3 + beta( uh( quadrature[i][0] ) ) * psi2 * psi3 );
      localMatrix[2][0] += wj * ( lambdas * gradPsi3 * gradPsi1 + sigma * psi3 * psi1 + beta( uh( quadrature[i][0] ) ) * psi3 * psi1 );
      localMatrix[2][1] += wj * ( lambdas * gradPsi3 * gradPsi2 + sigma * psi3 * psi2 + beta( uh( quadrature[i][0] ) ) * psi3 * psi2 );
      localMatrix[2][2] += wj * ( lambdas * gradPsi3 * gradPsi3 + sigma * psi3 * psi3 + beta( uh( quadrature[i][0] ) ) * psi3 * psi3 );
   }
}

void Element::BuildLocalB( )
{
   double **quadrature = basis.Quadrature( );
   double jacobian = coords[1] - coords[0];
   for ( int i = 0; i < 3; i++ )
   {
      double wj = quadrature[i][1] * jacobian;
      double psi1 = basis.basis1( quadrature[i][0] );
      double psi2 = basis.basis2( quadrature[i][0] );
      double psi3 = basis.basis3( quadrature[i][0] );
      double gradPsi1 = basis.basis1Grad( quadrature[i][0] );
      double gradPsi2 = basis.basis2Grad( quadrature[i][0] );
      double gradPsi3 = basis.basis3Grad( quadrature[i][0] );
      double lambdas = lambda( uh( quadrature[i][0] ) );
      double fVal = f[0] * psi1 + f[1] * psi2 + f[2] * psi3;
      localB[0] += wj * ( fVal * psi1 + theta( uh( quadrature[i][0] ) ) * psi1 + beta( uh( quadrature[i][0] ) ) * uBeta( uh( quadrature[i][0] ) ) * psi1 );
      localB[1] += wj * ( fVal * psi2 + theta( uh( quadrature[i][0] ) ) * psi2 + beta( uh( quadrature[i][0] ) ) * uBeta( uh( quadrature[i][0] ) ) * psi2 );
      localB[2] += wj * ( fVal * psi3 + theta( uh( quadrature[i][0] ) ) * psi3 + beta( uh( quadrature[i][0] ) ) * uBeta( uh( quadrature[i][0] ) ) * psi3 );
   }
}
