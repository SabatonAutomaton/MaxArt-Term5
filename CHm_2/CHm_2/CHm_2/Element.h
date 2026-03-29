#pragma once
#include "Basis.h"
#include <vector>
#include "Func.h"

class Element
{
public:
   Functions localFunctions;
   Functions* functions = &localFunctions;

   Element( ) {}
   Element( std::vector<double> coords_, int cond_, std::vector<double> q_, std::vector<double> _f, std::vector<double> _lambdaNodes, Functions* functions_ = nullptr )
   {
      coords = coords_;
      cond = cond_;
      localB.resize( 3 );
      f = _f;
      lambdaNodes = _lambdaNodes;
      q = q_;
      SetFunctions( functions_ );
   }

   void SetFunctions( Functions* functions_ );

   //параметры задачи
   double lambdaApprox( double x );
   double materialLambda( double u );
   double materialDLambda(double u);
   std::vector<double> lambdaNodes;
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

