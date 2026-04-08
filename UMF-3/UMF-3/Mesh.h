#pragma once
#include <vector>
class Mesh
{
   public:
      std::vector<double> meshX;
      int n_x;

      bool uniformX = true;

      void InputMeshes( );

      void CalculateUniformMesh( std::vector<double> &mesh, double a, double b, int n );
      void CalculateNonUniformMesh( std::vector<double> &mesh, double a, double b, int n, double q );
};

