#pragma once
#include "SLE.h"
#include "SLEAssistant.h"
#include <iostream>

class BCGSTAB
{
public:
   static std::vector<double> SolutionWithDiagonalConditioning( SLE& sle );
};
