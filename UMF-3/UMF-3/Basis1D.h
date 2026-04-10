#pragma once
class Basis1D
{
public:
   double Basis1( double x) const;
   double Basis2( double x ) const;
   double GradBasis1( double x ) const;
   double GradBasis2( double x ) const;
};

