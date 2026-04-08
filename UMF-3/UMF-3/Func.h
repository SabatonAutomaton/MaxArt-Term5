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
};

