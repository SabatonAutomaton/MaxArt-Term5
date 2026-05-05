#include "MatrixOperations.h"


Matrix MatrixOperations::MultMatScal( Matrix mat, double scal )
{
   Matrix res = mat;
   for ( int i = 0; i < mat.n; i++ )
   {
      res.di[i] *= scal;
      for ( int j = 0; j < 2; j++ )
      {
         res.ggl[j][i] *= scal;
         res.ggu[j][i] *= scal;
      }
   }

   return res;
}

std::vector<double> MatrixOperations::MultMatVec( Matrix mat, std::vector<double> vec )
{
   int n = mat.ggl[0].size( );
   std::vector<double> y( mat.n, 0.0 );

   for ( int i = 0; i < mat.n; i++ )
   {
      // диагональ
      y[i] += mat.di[i] * vec[i];

      int maxI = std::max( n - i, 0 );
      // левая часть (j < i)
      for ( int j = maxI; j < n; j++ )
      {
         y[j - n + i] += mat.ggu[i][j] * vec[i];
      }

      // правая часть (j > i)
      for ( int j = maxI; j < n; j++ )
      {
         y[i] += mat.ggl[i][j] * vec[j - n + i];
      }
   }

   return y;
}

std::vector<double> MatrixOperations::MultVecScal( std::vector<double> vec, double scal )
{
   std::vector<double> res = vec;
   for ( int i = 0; i < vec.size( ); i++ )
      res[i] *= scal;
   return res;
}

double MatrixOperations::DotProduct( std::vector<double> vec1, std::vector<double> vec2 )
{
   double res = 0;
   for ( int i = 0; i < vec1.size( ); i++ )
      res += vec1[i] * vec2[i];
   return res;
}

std::vector<double> MatrixOperations::AddVec( std::vector<double> vec1, std::vector<double> vec2 )
{
   std::vector<double> res = vec1;
   for ( int i = 0; i < vec1.size( ); i++ )
      res[i] += vec2[i];
   return res;
}
