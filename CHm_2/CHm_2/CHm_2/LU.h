#pragma once
#include "Matrix.h"
#include <vector>
class LU
{
public:
   LU( ) {}
   LU( Matrix _matrix, std::vector<double> _b )
   {
      matrix = _matrix;
      b = _b;
      y = std::vector<double>( matrix.n, 0.0 );
      q = std::vector<double>( matrix.n, 0.0 );
   }
   Matrix matrix;
   std::vector<double> b;
   std::vector<double> y;
   std::vector<double> q;

   void calcLU( );

   void calcY( );

   void calcQ( );

};

