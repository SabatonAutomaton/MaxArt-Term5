#include "LOS.h"

#include <cmath>
#include <iostream>

std::vector<double> LOS::SolutionWithoutConditioning(SLE sle)
{
    const int n = sle.matrix.n;

    const double fNorm =
        std::sqrt(
            SLEAssistant::DotProduct(
                sle.f,
                sle.f));

    if (fNorm < 1e-30)
    {
        std::cout << "Residual: 0 k: 0\n";
        return sle.x;
    }

    std::vector<double> r(n);
    std::vector<double> z(n);
    std::vector<double> p(n);
    std::vector<double> q(n);

    SLEAssistant::MatrixVectorProduct(
        sle.matrix,
        sle.x,
        q);

    for (int i = 0; i < n; ++i)
    {
        r[i] = sle.f[i] - q[i];
        z[i] = r[i];
    }

    SLEAssistant::MatrixVectorProduct(
        sle.matrix,
        z,
        p);

    int k = 0;
    double residual = 0.0;

    for (; k < sle.maxIter; ++k)
    {
        const double pp =
            SLEAssistant::DotProduct(p, p);

        if (std::abs(pp) < 1e-30)
            break;

        const double alpha =
            SLEAssistant::DotProduct(p, r) / pp;

        for (int i = 0; i < n; ++i)
        {
            sle.x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
        }

        residual =
            std::sqrt(
                SLEAssistant::DotProduct(r, r));

        if (residual / fNorm < sle.eps)
            break;

        SLEAssistant::MatrixVectorProduct(
            sle.matrix,
            r,
            q);

        const double beta =
            -SLEAssistant::DotProduct(p, q) / pp;

        for (int i = 0; i < n; ++i)
        {
            z[i] = r[i] + beta * z[i];
            p[i] = q[i] + beta * p[i];
        }
    }

    std::cout
        << "Residual: "
        << residual / fNorm
        << " k: "
        << k
        << "\n";

    return sle.x;
}

std::vector<double> LOS::SolutionWithDiagonalConditioning(SLE sle)
{
    const int n = sle.matrix.n;

    const double fNorm =
        std::sqrt(
            SLEAssistant::DotProduct(
                sle.f,
                sle.f));

    if (fNorm < 1e-30)
    {
        std::cout << "Residual: 0 k: 0\n";
        return sle.x;
    }

    std::vector<double> r(n);
    std::vector<double> z(n);
    std::vector<double> p(n);
    std::vector<double> q(n);
    std::vector<double> tmp(n);

    SLEAssistant::MatrixVectorProduct(
        sle.matrix,
        sle.x,
        tmp);

    for (int i = 0; i < n; ++i)
        r[i] = sle.f[i] - tmp[i];

    SLEAssistant::DiagonalPreconditionLOS(
        sle.matrix,
        r,
        z);

    SLEAssistant::MatrixVectorProduct(
        sle.matrix,
        z,
        p);

    double residual =
        std::sqrt(
            SLEAssistant::DotProduct(
                r,
                r));

    int k = 0;

    for (; k < sle.maxIter; ++k)
    {
        const double pp =
            SLEAssistant::DotProduct(
                p,
                p);

        if (std::abs(pp) < 1e-30)
            break;

        const double alpha =
            SLEAssistant::DotProduct(
                p,
                r) / pp;

        for (int i = 0; i < n; ++i)
        {
            sle.x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
        }

        residual =
            std::sqrt(
                SLEAssistant::DotProduct(
                    r,
                    r));

        if (residual / fNorm < sle.eps)
            break;

        SLEAssistant::DiagonalPreconditionLOS(
            sle.matrix,
            r,
            tmp);

        SLEAssistant::MatrixVectorProduct(
            sle.matrix,
            tmp,
            q);

        const double beta =
            -SLEAssistant::DotProduct(
                p,
                q) / pp;

        for (int i = 0; i < n; ++i)
        {
            z[i] = tmp[i] + beta * z[i];
            p[i] = q[i] + beta * p[i];
        }
    }

    std::cout
        << "Residual: "
        << residual / fNorm
        << " k: "
        << k
        << "\n";

    return sle.x;
}