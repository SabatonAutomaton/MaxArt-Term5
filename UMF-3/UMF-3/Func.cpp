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

   double lambdaCoeff( double z ) { ( void )z; return 2.0; }
   double dlambdaDz( ) { return 0.0; }
   double omegaCoeff( double z ) { ( void )z; return 5.0; }
   double sigmaCoeff( double z ) { ( void )z; return 0.7; }
   double xiCoeff( double z ) { ( void )z; return 0.4; }

   struct Linear3D
   {
      double c0;
      double cx;
      double cy;
      double cz;
   };

   Linear3D us3DCoefs( int variant )
   {

      return Linear3D{ -0.15, 0.35, 0.45, -0.25 };

   }

   Linear3D uc3DCoefs( int variant )
   {

      return Linear3D{ 0.3, -0.2, 0.15, 0.5 };

   }

   double evalLinear3D( const Linear3D &c, double x, double y, double z )
   {
      return c.c0 + c.cx * x + c.cy * y + c.cz * z;
   }

   double evalNonLinear3D( const Linear3D &c, double x, double y, double z )
   {
      return c.cx * x*x + c.cy * y*y + c.cz * z*z;
   }
}

void Func::SetDomain( double left, double right )
{
   gLeftZ = left;
   gRightZ = right;
}

double Func::u( double x )
{
   ( void )x;
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
   ( void )x;
   return 0.0;
}

double Func::uBetaS( double x )
{
   ( void )x;
   return 0.0;
}

double Func::uBetaC( double x )
{
   ( void )x;
   return 0.0;
}

double Func::boundaryValue( int number, double x )
{
   ( void )x;
   if ( number == 0 ) // Синусная компонента A^s
      return usExact( gLeftZ );
   if ( number == 1 ) // Косинусная компонента A^c
      return ucExact( gLeftZ );
   return 0.0;
}

double Func::lambda3D( int variant, double x, double y, double z ) const
{
   ( void )variant;
   ( void )x;
   ( void )y;
   ( void )z;
   return lambdaVal;
}

double Func::w3D( int variant, double x, double y, double z ) const
{
   ( void )variant;
   ( void )x;
   ( void )y;
   ( void )z;
   return wVal;
}

double Func::sigma3D( int variant, double x, double y, double z ) const
{
   ( void )variant;
   ( void )x;
   ( void )y;
   ( void )z;
   return sigmaVal;
}

double Func::xi3D( int variant, double x, double y, double z ) const
{
   ( void )variant;
   ( void )x;
   ( void )y;
   ( void )z;
   return xiVal;
}

double evalPoly3D(int p, double x, double y, double z)
{
    return std::pow(x, p)
        + std::pow(y, p)
        + std::pow(z, p);
}

double lapPoly3D(int p, double x, double y, double z)
{
    if (p < 2)
        return 0.0;

    return p * (p - 1) *
        (std::pow(x, p - 2)
            + std::pow(y, p - 2)
            + std::pow(z, p - 2));
}


double Func::us3DExact(
    int variant,
    double x,
    double y,
    double z) const
{
    if (variant == 2)
        return evalLinear3D(us3DCoefs(variant), x, y, z);

    if (variant >= 10)
    {
        const int p = variant - 9;
        return evalPoly3D(p, x, y, z);
    }

    return evalNonLinear3D(
        us3DCoefs(variant),
        x,
        y,
        z);
}

double Func::uc3DExact(
    int variant,
    double x,
    double y,
    double z) const
{
    if (variant == 2)
        return evalLinear3D(uc3DCoefs(variant), x, y, z);

    if (variant >= 10)
    {
        const int p = variant - 9;
        return 0.5 * evalPoly3D(p, x, y, z);
    }

    return evalNonLinear3D(
        uc3DCoefs(variant),
        x,
        y,
        z);
}

double Func::LapUs(
    int variant,
    double x,
    double y,
    double z) const
{
    if (variant == 2)
        return 0.0;

    if (variant >= 10)
    {
        const int p = variant - 9;
        return lapPoly3D(p, x, y, z);
    }

    return us3DCoefs(variant).cx
        + us3DCoefs(variant).cy
        + us3DCoefs(variant).cz;
}

double Func::LapUc(
    int variant,
    double x,
    double y,
    double z) const
{
    if (variant == 2)
        return 0.0;

    if (variant >= 10)
    {
        const int p = variant - 9;
        return 0.5 * lapPoly3D(p, x, y, z);
    }

    return uc3DCoefs(variant).cx
        + uc3DCoefs(variant).cy
        + uc3DCoefs(variant).cz;
}

double Func::fs3D( int variant, double x, double y, double z ) const
{
   const double lambdaVal = lambda3D( variant, x, y, z );
   const double wVal = w3D( variant, x, y, z );
   const double sigmaVal = sigma3D( variant, x, y, z );
   const double xiVal = xi3D( variant, x, y, z );
   const double lapUs = LapUs( variant, x, y, z );
   return -lambdaVal * lapUs - ( wVal * wVal ) * xiVal * us3DExact( variant, x, y, z ) - wVal * sigmaVal * uc3DExact( variant, x, y, z );
}

double Func::fc3D( int variant, double x, double y, double z ) const
{
   const double lambdaVal = lambda3D( variant, x, y, z );
   const double wVal = w3D( variant, x, y, z );
   const double sigmaVal = sigma3D( variant, x, y, z );
   const double xiVal = xi3D( variant, x, y, z );
   const double lapUc = LapUc( variant, x, y, z );
   return -lambdaVal * lapUc - ( wVal * wVal ) * xiVal * uc3DExact( variant, x, y, z ) + wVal * sigmaVal * us3DExact( variant, x, y, z );
}
