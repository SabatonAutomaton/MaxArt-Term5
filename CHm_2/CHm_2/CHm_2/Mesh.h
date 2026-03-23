#pragma once
#include <vector>
class Mesh
{
public:
   std::vector<double> meshX;
   std::vector<double> meshT;
   int n_x, n_t;

   bool uniformX = true;
   bool uniformT = true;

   void InputMeshes( );

   void CalculateUniformMesh( std::vector<double> &mesh, double a, double b, int n);
   void CalculateNonUniformMesh( std::vector<double> &mesh, double a, double b, int n, double q );
};

