#pragma once

#include "SparseMatrix.h"
#include <fstream>

class SLE
{
public:
   int                 maxIter = 0;
   double              eps = 0;
   int                 iterCount = 0;
   SparseMatrix        matrix;
   std::vector<double> f;
   std::vector<double> x;

   void Input( );

private:
   void InputProperty( );
   void InputSLE( );
};