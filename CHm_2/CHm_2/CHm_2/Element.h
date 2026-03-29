#pragma once
#include "Basis.h"
#include <vector>
#include "Func.h"

class Element
{
public:
   Functions functions;

   Element( ) {}
   Element( std::vector<double> coords_, int cond_, std::vector<double> q_, std::vector<double> _f, std::vector<double> _lambdaNodes )
   {
      coords = coords_;
      cond = cond_;
      localB.resize( 3 );
      f = _f;
      lambdaNodes = _lambdaNodes;
      q = q_;
   }

   //параметры задачи
   double lambda( double u );
   double dlambda(double u);
   std::vector<double> lambdaNodes;
   double sigma = 1;
   std::vector<double> f;

   double theta( double u );
   double beta( double u );
   double dbeta(double u);
   double uBeta( double u );
   double uh( double x );

   std::vector<double> coords;
   int cond = 0;
   std::vector<double> q;
   std::vector<std::vector<double>> localMatrix = std::vector<std::vector<double>>( 3, std::vector<double>( 3, 0.0 ) );
   std::vector<double> localB;

   Basis basis;

   void BuildLocalMatrix( );
   void BuildLocalB( );
   void BuildLocalMatrixNewton(std::vector<double> q_0);
   void BuildLocalBNewton(std::vector<double> q_0);
};

