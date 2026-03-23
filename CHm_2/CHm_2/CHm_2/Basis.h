#pragma once
class Basis
{
public:
   double basis1( double x );
   double basis1Grad( double x );
   double basis2( double x );
   double basis2Grad( double x );
   double basis3( double x );
   double basis3Grad( double x );

   double **Quadrature( );
};

