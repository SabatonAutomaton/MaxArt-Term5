#pragma once
#include "SparseMatrix.h"

// Неполное LU-разложение ILU(0):
// Fill-in за пределами профиля исходной матрицы отбрасывается.
// L — нижний треугольник в ggl (диагональ L = 1, не хранится)
// U — верхний треугольник в ggu
// D — диагональ U в di
class ILU : public SparseMatrix
{
public:
   std::vector<double> ggu;
   explicit ILU( const SparseMatrix &matrix );
};