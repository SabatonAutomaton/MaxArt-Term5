#include "Mesh.h"
#include "SimpleIteration.h"
#include <iostream>
#include "LU.h"

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


   SimpleIteration simpleIteration;
   simpleIteration.Input( );
   simpleIteration.Solve( );
   for ( int i = 0; i < simpleIteration.n; i++ )
   {
      std::cout << simpleIteration.q[i] << "\t";
   }
   int n = simpleIteration.n;
   //std::cout << std::endl;
   //std::vector<double> res = simpleIteration.ops.MultMatVec(simpleIteration.matrix,simpleIteration.q);
   //for ( int i = 0; i < n; i++ )
   //{
   //   std::cout << res[i] << "\t";
   //}
   //std::cout << std::endl;
   //for ( int i = 0; i < n; i++ )
   //{
   //   std::cout << simpleIteration.b[i] << "\t";
   //}
   //std::cout << std::endl;

   /*simpleIteration.BuildMatrix( );*/

   std::cout << "di = ";
   for ( int i = 0; i < n; i++ )
   {
      std::cout << simpleIteration.matrix.di[i] << "\t";
   }
   std::cout << std::endl;
   std::cout << "al = ";
   for ( int i = 0; i < n; i++ )
   {
      std::cout << simpleIteration.matrix.ggl[i][0] << "\t" << simpleIteration.matrix.ggl[i][1] << std::endl;
      std::cout << "\n     ";
   }
   std::cout << std::endl;
   std::cout << "au = ";
   for ( int i = 0; i < n; i++ )
   {
      std::cout << simpleIteration.matrix.ggu[i][0] << "\t" << simpleIteration.matrix.ggu[i][1] << std::endl;
      std::cout << "\n     ";
   }
   std::cout << std::endl << std::endl;

   std::cout << "Local matrixes:\n";
   for ( int k = 0; k < simpleIteration.elemCount; k++ )
   {
      std::cout << k << " element" << std::endl;
      for ( int i = 0; i < 3; i++ )
      {
         for ( int j = 0; j < 3; j++ )
            std::cout << simpleIteration.elements[k].localMatrix[i][j] << "\t";
         std::cout << std::endl;
      }
      std::cout << std::endl;
   }

   std::cout << "Local b:\n";
   for ( int k = 0; k < simpleIteration.elemCount; k++ )
   {
      std::cout << k << " element" << std::endl;
      for ( int i = 0; i < 3; i++ )
      {
         std::cout << simpleIteration.elements[k].localB[i] << "\t";
      }
      std::cout << std::endl;
   }

   //for ( int i = 0; i < n; i++ )
   //{
   //   for ( int j = 0; j < n; j++ )
   //   {
   //      if ( i == j )
   //         std::cout << simpleIteration.matrix.di[i] << "\t";
   //      else if ( j == i - 1 )  // вторая поддиагональ
   //         std::cout << simpleIteration.matrix.ggl[i][1] << "\t";
   //      else if ( j == i - 2 )  // первая поддиагональ
   //         std::cout << simpleIteration.matrix.ggl[i][0] << "\t";
   //      else if ( j == i + 2 && i != n - 1 )  // первая наддиагональ
   //         std::cout << simpleIteration.matrix.ggu[i + 1][0] << "\t";
   //      else if ( j == i + 1 && i != n - 1 )  // вторая наддиагональ
   //         std::cout << simpleIteration.matrix.ggu[i + 1][1] << "\t";
   //      else
   //         std::cout << "0\t";
   //   }
   //   std::cout << std::endl;
   //}

   std::cout << "\n\n";
   for ( int i = 0; i < n; i++ )
   {
      std::cout << simpleIteration.b[i] << "\t";
   }
   std::cout << "\n\n";

   //Matrix matrix;
   //matrix.n = 5;
   //matrix.di = { 10, 10, 10, 10, 10 };
   //matrix.ggl = { {0,0,0,0},{0,0,0,1},{0,0,2,3},{0,4,5,6},{7,8,9,11} };
   //matrix.ggu = { {0,0,0,0},{0,0,0,1},{0,0,2,3},{0,4,5,6},{7,8,9,11} };
   ///*matrix.ggl = { {1,2,4,7},{3,5,8,0},{6,9,0,0},{11,0,0,0},{0,0,0,0} };
   //matrix.ggu = { {1,2,4,7},{3,5,8,0},{6,9,0,0},{11,0,0,0},{0,0,0,0} };*/
   //std::vector<double> b = { 24,27,30,36,45 };
   //LU lu( matrix, b );
   //lu.calcLU( );
   //lu.calcY( );
   //lu.calcQ( );
   //for ( int i = 0; i < matrix.n; i++ )
   //{
   //   std::cout << lu.q[i] << "\t";
   //}

   //Matrix A;
   //A.n = 5;
   //A.di = { 10,10,10,10,10 };
   //A.ggl = { {0,0},
   //          {0,1},
   //          {2,3},
   //          {5,6},
   //          {9,11} };
   //A.ggu = A.ggl;
   //A.ggu[1][1] = 0;
   //A.ggu[2][0] = 0;
   ///*A.ggl[1][1] = 0;
   //A.ggl[1][0] = 0;*/

   //std::vector<double> x = { 1,1,1,1,1 };

   //auto y = simpleIteration.ops.MultMatVec( A, x );

   //for ( int i = 0; i < 5; i++ )
   //{
   //   std::cout << y[i] << "\t";
   //}
}