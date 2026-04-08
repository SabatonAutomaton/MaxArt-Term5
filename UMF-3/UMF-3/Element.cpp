#include "Element.h"

void Element::CalcCond( )
{
   for ( int i = 0; i < 2; i++ )
   {
      localB[2 * i] += func.thetaS(coords[i] ) + uBetaS;
      localB[2 * i + 1] += func.thetaC( coords[i] ) + uBetaC;
      for ( int j = 0; j < 2; j++ )
      {
         localMatrix[2 * i][2 * j] += beta;
         localMatrix[2 * i + 1][2 * j + 1] += beta;
      }
   }
}

double Element::Approx( std::vector<double> val, double x )
{
   return val[0] * basis.Basis1( x ) + val[1] * basis.Basis2( x );
}

void Element::BuildLocalMatrix( )
{
   std::vector<double> psi( 2, 0.0 );
   std::vector<double> gradPsi( 2, 0.0 );
   std::vector<std::vector<double>> quadratures = quad.gauss3( );
   const double h = coords[1] - coords[0];

   for ( auto q : quadratures )
   {
      const double wj = q[1] * h;

      psi[0] = basis.Basis1( q[0] );
      psi[1] = basis.Basis2( q[0] );
      // Производная в физической координате: dphi/dx = dphi/dxi * dxi/dx = dphi/dxi / h
      gradPsi[0] = basis.GradBasis1( q[0] ) / h;
      gradPsi[1] = basis.GradBasis2( q[0] ) / h;
      double lambda = Approx( lambdaNodes, q[0] );
      double w = Approx( wNodes, q[0] );
      double sigma = Approx( sigmaNodes, q[0] );
      double xi = Approx( xiNodes, q[0] );

      for ( int i = 0; i < 2; i++ )
      {
         for ( int j = 0; j < 2; j++ )
         {
            const double a = wj * ( lambda * gradPsi[i] * gradPsi[j] - w * w * xi * psi[i] * psi[j] );
            const double b = wj * w * sigma * psi[i] * psi[j];

            localMatrix[2 * i][2 * j] += a;
            localMatrix[2 * i][2 * j + 1] += -b;
            localMatrix[2 * i + 1][2 * j] += b;
            localMatrix[2 * i + 1][2 * j + 1] += a;
         }
      }

   }
}

void Element::BuildLocalB( )
{
   std::vector<double> psi( 2, 0.0 );
   std::vector<double> gradPsi( 2, 0.0 );
   std::vector<std::vector<double>> quadratures = quad.gauss3( );
   double jacobian = coords[1] - coords[0];

   for ( auto q : quadratures )
   {
      double wj = q[1] * jacobian;

      psi[0] = basis.Basis1( q[0] );
      psi[1] = basis.Basis2( q[0] );
      gradPsi[0] = basis.GradBasis1( q[0] );
      gradPsi[1] = basis.GradBasis2( q[0] );
      double fs = Approx( fsNodes, q[0] );
      double fc = Approx( fcNodes, q[0] );

      for ( int i = 0; i < 2; i++ )
      {
         localB[2 * i] += wj * fs * psi[i];
         localB[2 * i + 1] += wj * fc * psi[i];
      }
   }
}
