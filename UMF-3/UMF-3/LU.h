#pragma once

#include "SparseMatrix.h"

class LU : public SparseMatrix
{
public:
   std::vector<double> ggu;

   LU( const SparseMatrix &matrix );
};