#include "Mesh.h"
#include <iostream>

int main( )
{
   Mesh mesh;
   mesh.InputMeshes();
   for ( int i = 0; i < mesh.n_x; i++ )
   {
      std::cout << mesh.meshX[i]<<std::endl;
   }
   for ( int i = 0; i < mesh.n_t; i++ )
   {
      std::cout << mesh.meshT[i] << std::endl;
   }
   return 0;
}