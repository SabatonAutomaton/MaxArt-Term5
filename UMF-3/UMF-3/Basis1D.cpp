#include "Basis1D.h"

double Basis1D::Basis1( double x )
{
   return 1-x;
}

double Basis1D::Basis2( double x )
{
   return x;
}

double Basis1D::GradBasis1( double x )
{
   return -1;
}

double Basis1D::GradBasis2( double x )
{
   return 1;
}
