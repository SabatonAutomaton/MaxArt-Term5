#!/bin/bash

# Группа 1: λ=100, ω=100, σ=1
echo "=== Группа 1: λ=100, ω=100, σ=1 ==="
for chi in 1e-11 1000 10000 100000 800000; do
    echo "100 100 1 $chi" > coefficients.txt

    echo "MSG" > solver.txt
    echo "Тест: χ=$chi, MSG"
    ./x64/Release/UMF-3.exe 2 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"

    echo "LU" > solver.txt
    echo "Тест: χ=$chi, LU"
    ./x64/Release/UMF-3.exe 2 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"
    echo ""
done

# Группа 2: λ=1000, ω=0.0001, σ=1
echo "=== Группа 2: λ=1000, ω=0.0001, σ=1 ==="
for chi in 1e-11 0.01 1 100 1e4 1e6 1e9; do
    echo "1000 0.0001 1 $chi" > coefficients.txt

    echo "MSG" > solver.txt
    echo "Тест: χ=$chi, MSG"
    ./x64/Release/UMF-3.exe 2 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"

    echo "LU" > solver.txt
    echo "Тест: χ=$chi, LU"
    ./x64/Release/UMF-3.exe 2 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"
    echo ""
done

# Группа 3: λ=1000, ω=100, σ=0
echo "=== Группа 3: λ=1000, ω=100, σ=0 ==="
for chi in 1e-11 100 10000 1e6 1e8; do
    echo "1000 100 0 $chi" > coefficients.txt

    echo "MSG" > solver.txt
    echo "Тест: χ=$chi, MSG"
    ./x64/Release/UMF-3.exe 2 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"

    echo "LU" > solver.txt
    echo "Тест: χ=$chi, LU"
    ./x64/Release/UMF-3.exe 2 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"
    echo ""
done

# Группа 4: λ=1000, ω=100, σ=1
echo "=== Группа 4: λ=1000, ω=100, σ=1 ==="
for chi in 8.81e-12 1e-12 1e-11 1e-10; do
    echo "1000 100 1 $chi" > coefficients.txt

    echo "MSG" > solver.txt
    echo "Тест: χ=$chi, MSG"
    ./x64/Release/UMF-3.exe 2 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"

    echo "LU" > solver.txt
    echo "Тест: χ=$chi, LU"
    ./x64/Release/UMF-3.exe 2 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"
    echo ""
done
