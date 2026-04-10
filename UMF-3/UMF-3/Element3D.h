#pragma once
#include "Basis3D.h"
#include "Func.h"
#include "Quadratures.h"
#include <array>
#include <vector>

class Element3D
{
public:
   using Point3 = std::array<double, 3>;

   Element3D( ) = default;
   Element3D(
      const std::array<Point3, 8> &_coords,
      const std::vector<double> &_fsNodes,
      const std::vector<double> &_fcNodes,
      const std::vector<double> &_lambdaNodes,
      const std::vector<double> &_wNodes,
      const std::vector<double> &_sigmaNodes,
      const std::vector<double> &_xiNodes )
      : coords( _coords ),
      fsNodes( _fsNodes ),
      fcNodes( _fcNodes ),
      lambdaNodes( _lambdaNodes ),
      wNodes( _wNodes ),
      sigmaNodes( _sigmaNodes ),
      xiNodes( _xiNodes )
   {
      localMatrix = std::vector<std::vector<double>>( 16, std::vector<double>( 16, 0.0 ) );
      localB = std::vector<double>( 16, 0.0 );
   }

   std::array<Point3, 8> coords {};
   std::vector<std::vector<double>> localMatrix;
   std::vector<double> localB;
   std::vector<double> fsNodes;
   std::vector<double> fcNodes;
   std::vector<double> lambdaNodes;
   std::vector<double> wNodes;
   std::vector<double> sigmaNodes;
   std::vector<double> xiNodes;

   double Approx( const std::vector<double> &val, const std::array<double, 3> &rst ) const;
   void BuildLocalMatrix( );
   void BuildLocalB( );

private:
   Func func;
   Basis3D basis;
   Quadratures quad;

   static double Det3( const std::array<std::array<double, 3>, 3> &m );
   static std::array<std::array<double, 3>, 3> Inverse3( const std::array<std::array<double, 3>, 3> &m );
};
