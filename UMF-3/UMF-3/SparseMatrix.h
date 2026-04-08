#pragma once
#include <vector>
class SparseMatrix
{
public:
   SparseMatrix( ) {}
   SparseMatrix( int _n )
   {
      n = _n;
      di = std::vector( n, 0.0 );
      ig = std::vector<int>( n + 1, 0 );
      jg.clear( );
      ggl.clear( );
      ggu.clear( );
   }

   int n;

   std::vector<double> di;
   std::vector<int> ig;
   std::vector<int> jg;
   std::vector<double> ggl;
   std::vector<double> ggu;
};

