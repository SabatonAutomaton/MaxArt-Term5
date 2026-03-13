#pragma once
class Mesh
{
public:
   double *meshX = nullptr;
   double *meshT = nullptr;
   int n_x, n_t;

   bool uniformX = true;
   bool uniformT = true;

   void InputMeshes( );

   void CalculateUniformMesh( double *mesh, double a, double b, int n);
   void CalculateNonUniformMesh( double *mesh, double a, double b, int n, double q );
};

