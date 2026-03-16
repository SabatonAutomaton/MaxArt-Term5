#include "MatrixOperations.h"

double **MatrixOperations::AddMatrix( double **mat1, double **mat2 )
{
   double **res = new double *[n];
   for ( int i = 0; i < n; i++ )
   {
      res[i] = new double[n];
      for ( int j = 0; j < n; j++ )
         res[i][j] = mat1[i][j] + mat2[i][j];
   }
   return res;
}

double **MatrixOperations::MultVecVecT( double *vec1, double *vec2 )
{
   double **res = new double *[n];
   for ( int i = 0; i < n; i++ )
   {
      res[i] = new double[n];
      for ( int j = 0; j < n; j++ )
         res[i][j] = vec1[i] * vec2[j];
   }
   return res;
}

double **MatrixOperations::MultMatScal( double **mat, double scal )
{
   double **res = new double *[n];
   for ( int i = 0; i < n; i++ )
   {
      res[i] = new double[n];
      for ( int j = 0; j < n; j++ )
         res[i][j] = mat[i][j] * scal;
   }
   return res;
}

double *MatrixOperations::MultMatVec( double **mat, double *vec )
{
   double *res = new double[n] {0};
   for ( int i = 0; i < n; i++ )
      for ( int j = 0; j < n; j++ ) 
         res[i] += mat[i][j] * vec[j];
   return res;
}
