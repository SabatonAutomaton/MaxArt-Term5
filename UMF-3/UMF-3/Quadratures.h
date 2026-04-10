#pragma once
#include <vector>
class Quadratures
{
public:
   std::vector<std::vector<double>> gauss3();
   std::vector<std::vector<double>> gaussHex2();
};

