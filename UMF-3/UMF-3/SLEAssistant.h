#pragma once
#include "SLE.h"
#include "LU.h"
#include "LLt.h"

class SLEAssistant
{
public:
   static double DotProduct(
      const std::vector<double> &v1,
      const std::vector<double> &v2
   );
   static void MatrixVectorProduct(
      const SparseMatrix &m,
      const std::vector<double> &v,
      std::vector<double> &res
   );
   static void LUDirect(
      const LU &lu,
      const std::vector<double> &v,
      std::vector<double> &res
   );
   static void LUReverse(
      const LU &lu,
      const std::vector<double> &v,
      std::vector<double> &res
   );
   static void DiagonalPreconditionMSG(
      SparseMatrix &M,
      std::vector<double> &v,
      std::vector<double> &res
   );
   static void DiagonalPreconditionLOS(
      SparseMatrix &M,
      std::vector<double> &v,
      std::vector<double> &res
   );
   static void LLtDirect(
      const LLt &llt,
      const std::vector<double> &v,
      std::vector<double> &res
   );
   static void LLtReverse(
      const LLt &llt,
      const std::vector<double> &v,
      std::vector<double> &res
   );
};