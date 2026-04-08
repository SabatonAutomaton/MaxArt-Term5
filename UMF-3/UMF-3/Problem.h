#pragma once
#include "Element.h"
#include "SparseMatrix.h"
#include "Mesh.h"
#include <string>
#include <vector>
class Problem
{
public:
   void Run( );
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
   SolveMethod solveMethod = SolveMethod::LU;

   int elemCount = 0;
   std::vector<Element> elements;

   void Input( );
   void BuildMatrix( );
   void BuildB( );
   void Solve( );
   void ApplyDirichlet( int dof, double value );
   void PrintResult( ) const;
   void PrintSystem( const std::vector<std::vector<double>> &A, const std::vector<double> &rhs ) const;

   Mesh mesh;
};

