#!/bin/bash

echo "3" > test.txt

# Группа 1: Варьируется λ (ω=100, σ=1, χ=1e-11)
echo "=== Группа 1: Варьируется λ (ω=100, σ=1, χ=1e-11) ==="
for lambda in 100 1000 10000 100000 800000; do
    echo "$lambda 100 1 1e-11" > coefficients.txt
    echo "BCGSTAB" > solver.txt
    echo "Тест: λ=$lambda, BCGSTAB"
    ./x64/Release/UMF-3.exe 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"
    echo ""
done

# Группа 2: Варьируется ω (λ=1000, σ=1, χ=1e-11)
echo "=== Группа 2: Варьируется ω (λ=1000, σ=1, χ=1e-11) ==="
for omega in 0.0001 0.01 1 100 10000 1000000 1000000000; do
    echo "1000 $omega 1 1e-11" > coefficients.txt
    echo "BCGSTAB" > solver.txt
    echo "Тест: ω=$omega, BCGSTAB"
    ./x64/Release/UMF-3.exe 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"
    echo ""
done

# Группа 3: Варьируется σ (λ=1000, ω=100, χ=1e-11)
echo "=== Группа 3: Варьируется σ (λ=1000, ω=100, χ=1e-11) ==="
for sigma in 0 100 10000 1000000 100000000; do
    echo "1000 100 $sigma 1e-11" > coefficients.txt
    echo "BCGSTAB" > solver.txt
    echo "Тест: σ=$sigma, BCGSTAB"
    ./x64/Release/UMF-3.exe 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"
    echo ""
done

# Группа 4: Варьируется χ (λ=1000, ω=100, σ=1)
echo "=== Группа 4: Варьируется χ (λ=1000, ω=100, σ=1) ==="
for chi in 8.81e-12 1e-12 1e-11 1e-10; do
    echo "1000 100 1 $chi" > coefficients.txt
    echo "BCGSTAB" > solver.txt
    echo "Тест: χ=$chi, BCGSTAB"
    ./x64/Release/UMF-3.exe 2>&1 | grep -E "(Time:|u\* - u\|L2 =)"
    echo ""
done

echo "=== Все тесты BCGSTAB завершены ==="
