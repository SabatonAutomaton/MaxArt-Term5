#pragma once
#include <vector>
#include "Matrix.h"
class MatrixOperations
{
public:
   Matrix MultMatScal( Matrix mat, double scal );
   std::vector<double> MultMatVec( Matrix mat, std::vector<double> vec );
   std::vector<double> MultVecScal( std::vector<double> vec, double scal );
   double DotProduct( std::vector<double> vec1, std::vector<double> vec2 );
   std::vector<double> AddVec( std::vector<double> vec1, std::vector<double> vec2 );
};