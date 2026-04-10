#pragma once
#include "Mesh.h"
#include <array>
#include <vector>

class Mesh3D : public Mesh
{
public:
   struct Node3D
   {
      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
   };

   struct Hex8
   {
      std::array<int, 8> v {};
   };

   std::vector<double> meshY;
   std::vector<double> meshZ;
   int n_y = 0;
   int n_z = 0;
   bool uniformY = true;
   bool uniformZ = true;

   std::vector<Node3D> nodes;
   std::vector<Hex8> hexes;

   void InputMeshes( );
   void BuildStructuredHexMesh( );
};
