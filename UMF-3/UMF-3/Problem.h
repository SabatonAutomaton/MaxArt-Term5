#pragma once
#include "Element.h"
#include "Element3D.h"
#include "Mesh3D.h"
#include "SparseMatrix.h"
#include "Mesh.h"
#include <string>
#include <vector>
class Problem
{
public:
   enum class Dimension
   {
      OneD = 1,
      ThreeD = 3
   };

   void Run( int testCase );
   const std::vector<double> &GetSolution( ) const { return q; }

private:
   enum class SolveMethod
   {
      LU,
      MSGDiagonal
   };

   Func func;

   std::vector<double> q;

   int n = 0;
   SparseMatrix matrix;

   std::vector<double> b;
   std::vector<double> boundaryB;
   std::vector<char> fixedDofs;
   SolveMethod solveMethod = SolveMethod::LU;
   Dimension dimension = Dimension::OneD;

   int elemCount = 0;
   std::vector<Element> elements;
   std::vector<Element3D> elements3D;
   std::vector<double> q3d;
   std::vector<double> usExact3D;
   std::vector<double> ucExact3D;
   bool use3DLinearVariant2 = false;

   void Input( );
   void BuildMatrix( );
   void BuildB( );
   void Solve( );
   void ApplyDirichlet( int dof, double value );
   void PrintResult( ) const;

   void Input3D( bool quadraticExact );
   void BuildMatrix3D( );
   void BuildB3D( );
   void ApplyDirichlet3D( int dof, double value );
   void Solve3D( );
   void PrintResult3D( ) const;
   void Run3DHarmonicTest( bool quadraticExact );

   Mesh mesh;
   Mesh3D mesh3d;
};

