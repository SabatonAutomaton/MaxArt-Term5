#pragma once
class Func
{
public:
   void SetDomain( double left, double right );
   double u( double x );
   double fc( double x );
   double fs( double x );
   double lambda( double x );
   double w( double x );
   double sigma( double x );
   double xi( double x );
   double thetaS( double x );
   double thetaC( double x );
   double beta( double x );
   double uBetaS( double x );
   double uBetaC( double x );
   double boundaryValue( int number, double x );

   // 3D тестовые функции для гармонической системы.
   double lambdaVal;
   double wVal;
   double sigmaVal;
   double xiVal;
   double lambda3D( int variant, double x, double y, double z ) const;
   double w3D( int variant, double x, double y, double z ) const;
   double sigma3D( int variant, double x, double y, double z ) const;
   double xi3D( int variant, double x, double y, double z ) const;
   double us3DExact( int variant, double x, double y, double z ) const;
   double uc3DExact( int variant, double x, double y, double z ) const;
   double LapUs( int variant, double x, double y, double z ) const;
   double LapUc( int variant, double x, double y, double z ) const;
   double fs3D( int variant, double x, double y, double z ) const;
   double fc3D( int variant, double x, double y, double z ) const;
};

