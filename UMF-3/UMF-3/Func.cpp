#include "Func.h"
#include <cmath>

namespace
{
   double gLeftZ = 0.0;
   double gRightZ = 1.0;

   double len( )
   {
      const double l = gRightZ - gLeftZ;
      return ( std::abs( l ) > 0.0 ) ? l : 1.0;
   }

   double tNorm( double z )
   {
      return ( z - gLeftZ ) / len( );
   }

   // Точное линейное решение теста:
   // us(z) = 0.2 + 0.6 * t, uc(z) = -0.4 + 0.8 * t, t в [0,1].
   double usExact( double z ) { return 0.2 + 0.6 * tNorm( z ); }
   double ucExact( double z ) { return -0.4 + 0.8 * tNorm( z ); }
   double dusExact( ) { return 0.6 / len( ); }
   double ducExact( ) { return 0.8 / len( ); }

   double lambdaCoeff( double z ) { (void)z; return 2.0; }
   double dlambdaDz( ) { return 0.0; }
   double omegaCoeff( double z ) { (void)z; return 5.0; }
   double sigmaCoeff( double z ) { (void)z; return 0.7; }
   double xiCoeff( double z ) { (void)z; return 0.4; }
}

void Func::SetDomain( double left, double right )
{
   gLeftZ = left;
   gRightZ = right;
}

double Func::u( double x )
{
    (void)x;
    return 0.0;
}

double Func::fc( double x )
{
   const double lPrime = dlambdaDz( );
   const double duc = ducExact( );
   const double w = omegaCoeff( x );
   const double sigma = sigmaCoeff( x );
   const double xi = xiCoeff( x );

   // -(lambda * uc')' - w^2 * xi * uc + w * sigma * us = fc
   return -lPrime * duc - ( w * w ) * xi * ucExact( x ) + w * sigma * usExact( x );
}

double Func::fs( double x )
{
   const double lPrime = dlambdaDz( );
   const double dus = dusExact( );
   const double w = omegaCoeff( x );
   const double sigma = sigmaCoeff( x );
   const double xi = xiCoeff( x );

   // -(lambda * us')' - w^2 * xi * us - w * sigma * uc = fs
   return -lPrime * dus - ( w * w ) * xi * usExact( x ) - w * sigma * ucExact( x );
}

double Func::lambda( double x )
{
   return lambdaCoeff( x );
}

double Func::w( double x )
{
   return omegaCoeff( x );
}

double Func::sigma( double x )
{
   return sigmaCoeff( x );
}

double Func::xi( double x )
{
   return xiCoeff( x );
}

double Func::thetaS( double x )
{
   return lambdaCoeff( x ) * dusExact( );
}

double Func::thetaC( double x )
{
   return lambdaCoeff( x ) * ducExact( );
}

double Func::beta( double x )
{
   (void)x;
   return 0.0;
}

double Func::uBetaS( double x )
{
   (void)x;
   return 0.0;
}

double Func::uBetaC( double x )
{
   (void)x;
   return 0.0;
}

double Func::boundaryValue( int number, double x )
{
   (void)x;
   if ( number == 0 ) // Синусная компонента A^s
      return usExact( gLeftZ );
   if ( number == 1 ) // Косинусная компонента A^c
      return ucExact( gLeftZ );
   return 0.0;
}

