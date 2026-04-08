#pragma once
#include "Basis1D.h"
#include "Func.h"
#include "Quadratures.h"
#include <vector>
class Element
{
public:
   Element( ) {}

   Element(std::vector<double> _coords, std::vector<double> fs, std::vector<double> fc )
   {
      coords = _coords;
      fsNodes = fs;
      fcNodes = fc;
      localMatrix = std::vector<std::vector<double>>( 4, std::vector<double>( 4, 0.0 ) );
      localB = std::vector<double>( 4, 0.0 );
       lambdaNodes = { func.lambda( coords[0] ), func.lambda( coords[1] ) };
      wNodes = { func.w( coords[0] ), func.w( coords[1] ) };
      sigmaNodes = { func.sigma( coords[0] ), func.sigma( coords[1] ) };
      xiNodes = { func.xi( coords[0] ), func.xi( coords[1] ) };
      beta = func.beta( 0 );

   }
   std::vector<double> lambdaNodes;
   std::vector<double> wNodes;
   std::vector<double> sigmaNodes;
   std::vector<double> xiNodes;

   Func func;
   double beta;
   double thetaS = 0.0;
   double thetaC = 0.0;
   double uBetaS = 0.0;
   double uBetaC = 0.0;
   void CalcCond( );

   double Approx( std::vector<double> val, double x );

   Quadratures quad;
   Basis1D basis;
   std::vector<std::vector<double>> localMatrix;
   std::vector<double> fsNodes;
   std::vector<double> fcNodes;
   std::vector<double> localB;
   std::vector<double> coords;

   void BuildLocalMatrix( );
   void BuildLocalB( );
};

