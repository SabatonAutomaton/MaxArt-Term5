#include "Basis.h"

double Basis::basis1( double x )
{
   return 2 * ( x - 0.5 ) * ( x - 1 );
}

double Basis::basis1Grad( double x )
{
   return 4 * x - 3;
}

double Basis::basis2( double x )
{
   return -4 * x * ( x - 1 );
}

double Basis::basis2Grad( double x )
{
   return -4 * ( 2 * x - 1 );
}

double Basis::basis3( double x )
{
   return 2 * x * ( x - 0.5 );
}

double Basis::basis3Grad( double x )
{
   return 4 * x - 1;
}
