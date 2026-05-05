#include <iostream>
#include <vector>
#include <cmath>
#include "BCGSTAB.h"
#include "SLE.h"
#include "SparseMatrix.h"

int main()
{
    // Создаем простую несимметричную матрицу 5x5
    // Решение: x = [1, 1, 1, 1, 1]
    const int n = 5;

    SparseMatrix A(n);

    // Заполняем матрицу
    // Диагональ
    A.di = {4.0, 5.0, 6.0, 5.0, 4.0};

    // Профиль матрицы (нижний треугольник)
    A.ig = {0, 0, 1, 3, 6, 10};

    // Индексы столбцов
    A.jg = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3};

    // Нижний треугольник
    A.ggl = {0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0};

    // Верхний треугольник (несимметричная часть)
    A.ggu = {0.0, -1.0, 0.0, -1.0, -1.0, 0.0, -1.0, -1.0, -1.0, 0.0};

    // Правая часть: b = A * [1, 1, 1, 1, 1]
    std::vector<double> b(n);

    // Вычисляем b = A * x, где x = [1, 1, 1, 1, 1]
    for (int i = 0; i < n; ++i)
    {
        b[i] = A.di[i]; // диагональный элемент

        // Нижний треугольник
        for (int k = A.ig[i]; k < A.ig[i + 1]; ++k)
        {
            b[i] += A.ggl[k];
        }

        // Верхний треугольник
        for (int j = i + 1; j < n; ++j)
        {
            for (int k = A.ig[j]; k < A.ig[j + 1]; ++k)
            {
                if (A.jg[k] == i)
                {
                    b[i] += A.ggu[k];
                    break;
                }
            }
        }
    }

    std::cout << "Тест BCGSTAB на малой матрице 5x5\n";
    std::cout << "Ожидаемое решение: x = [1, 1, 1, 1, 1]\n\n";

    std::cout << "Правая часть b:\n";
    for (int i = 0; i < n; ++i)
        std::cout << "b[" << i << "] = " << b[i] << "\n";
    std::cout << "\n";

    // Решаем систему BCGSTAB
    SLE sle;
    sle.matrix = A;
    sle.f = b;
    sle.x.assign(n, 0.0);
    sle.maxIter = 1000;
    sle.eps = 1e-12;

    std::cout << "Запуск BCGSTAB...\n";
    std::vector<double> x = BCGSTAB::SolutionWithDiagonalConditioning(sle);

    std::cout << "\nПолученное решение:\n";
    for (int i = 0; i < n; ++i)
        std::cout << "x[" << i << "] = " << x[i] << "\n";

    // Проверяем точность
    double error = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double diff = x[i] - 1.0;
        error += diff * diff;
    }
    error = std::sqrt(error);

    std::cout << "\nОшибка ||x - [1,1,1,1,1]|| = " << error << "\n";

    if (error < 1e-10)
        std::cout << "✓ ТЕСТ ПРОЙДЕН\n";
    else
        std::cout << "✗ ТЕСТ НЕ ПРОЙДЕН\n";

    return 0;
}
