#pragma once
class MatrixOperations
{
public:
   int n = 2;
   double **AddMatrix( double **mat1, double **mat2 );
   double **MultVecVecT(double *vec1, double *vec2 );
   double **MultMatScal(double **mat, double scal );
   double *MultMatVec( double **mat, double *vec );
};

