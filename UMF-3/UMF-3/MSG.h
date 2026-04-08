#pragma once
#include "SLE.h"
#include "SLEAssistant.h"
#include <iostream>

class MSG
{
public:
   static std::vector<double> SolutionWithoutConditioning( SLE sle );
   static std::vector<double> SolutionWithDiagonalConditioning( SLE sle );
   static std::vector<double> SolutionWithLUConditioning( SLE sle );
   static std::vector<double> SolutionWithLLtConditioning( SLE sle );
};
