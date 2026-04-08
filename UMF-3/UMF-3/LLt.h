#pragma once

#include "SparseMatrix.h"

class LLt : public SparseMatrix
{
public:
   LLt( const SparseMatrix &matrix );
};