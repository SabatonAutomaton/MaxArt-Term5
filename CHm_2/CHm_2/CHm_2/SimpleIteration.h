#pragma once
#include "Basis.h"

class SimpleIteration
{
public:
   //начальное приблежение
   double *q0 = nullptr;

   //матрица
   double **al = nullptr;
   double **au = nullptr;
   double *di = nullptr;

   //вектор правой части
   double *b = nullptr;

   void BuildMatrix( );
   void BuildB( );

   Basis basis;
};