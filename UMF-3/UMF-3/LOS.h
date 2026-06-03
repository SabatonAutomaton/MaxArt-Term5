#pragma once

#include "SLE.h"
#include "SLEAssistant.h"
#include <vector>

class LOS
{
public:
	static std::vector<double> SolutionWithoutConditioning(SLE& sle);
	static std::vector<double> SolutionWithDiagonalConditioning(SLE& sle);
};