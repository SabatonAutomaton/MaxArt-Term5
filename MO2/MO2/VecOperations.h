#pragma once
#include <iostream>
class VecOperations
{
public:
   int n= 2;
   double *MultVecScal( double *vec, double scal );
   double *AddVec( double *vec1, double *vec2 );
   double DotProduct( double *vec );
};

