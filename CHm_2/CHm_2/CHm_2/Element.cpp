#include "Element.h"

void Element::SetFunctions( Functions* functions_ )
{
   functions = functions_ != nullptr ? functions_ : &localFunctions;
}

double Element::lambdaApprox( double x )
{
   return lambdaNodes[0] * basis.basis1( x ) + lambdaNodes[1] * basis.basis2( x ) + lambdaNodes[2] * basis.basis3( x );
}

double Element::materialLambda( double u )
{
   return functions->lambda( u );
}

double Element::materialDLambda( double u )
{
   return functions->dlambda( u );
}

double Element::theta( double u )
{
   return functions->theta( u );
}

double Element::beta( double u )
{
   return functions->beta( u );
}

double Element::dbeta( double u )
{
   return functions->dbeta( u );
}

double Element::uBeta( double u )
{
   return functions->uBeta( u );
}

double Element::uh( double x )
{
   return q[0] * basis.basis1( x ) + q[1] * basis.basis2( x ) + q[2] * basis.basis3( x );
}

double Element::uhPrevTime( double x )
{
   return qPrevTime[0] * basis.basis1( x ) + qPrevTime[1] * basis.basis2( x ) + qPrevTime[2] * basis.basis3( x );
}

void Element::BuildLocalMatrix( )
{
   localMatrix = std::vector<std::vector<double>>( 3, std::vector<double>( 3, 0.0 ) );
   double **quadrature = basis.Quadrature( );
   double jacobian = coords[1] - coords[0];
   for ( int i = 0; i < 3; i++ )
   {
      double wj = quadrature[i][1] * jacobian;
      double psi1 = basis.basis1( quadrature[i][0] );
      double psi2 = basis.basis2( quadrature[i][0] );
      double psi3 = basis.basis3( quadrature[i][0] );
      double gradPsi1 = basis.basis1Grad( quadrature[i][0] ) / jacobian;
      double gradPsi2 = basis.basis2Grad( quadrature[i][0] ) / jacobian;
      double gradPsi3 = basis.basis3Grad( quadrature[i][0] ) / jacobian;
      double sigmaVal = useTimeDependent ? 0.0 : functions->sigma( uh( quadrature[i][0] ) );
      double lambdas = lambdaApprox( quadrature[i][0] );
      localMatrix[0][0] += wj * ( lambdas * gradPsi1 * gradPsi1 + sigmaVal * psi1 * psi1 + beta( uh( quadrature[i][0] ) ) * psi1 * psi1 );
      localMatrix[0][1] += wj * ( lambdas * gradPsi1 * gradPsi2 + sigmaVal * psi1 * psi2 + beta( uh( quadrature[i][0] ) ) * psi1 * psi2 );
      localMatrix[0][2] += wj * ( lambdas * gradPsi1 * gradPsi3 + sigmaVal * psi1 * psi3 + beta( uh( quadrature[i][0] ) ) * psi1 * psi3 );
      localMatrix[1][0] += wj * ( lambdas * gradPsi2 * gradPsi1 + sigmaVal * psi2 * psi1 + beta( uh( quadrature[i][0] ) ) * psi2 * psi1 );
      localMatrix[1][1] += wj * ( lambdas * gradPsi2 * gradPsi2 + sigmaVal * psi2 * psi2 + beta( uh( quadrature[i][0] ) ) * psi2 * psi2 );
      localMatrix[1][2] += wj * ( lambdas * gradPsi2 * gradPsi3 + sigmaVal * psi2 * psi3 + beta( uh( quadrature[i][0] ) ) * psi2 * psi3 );
      localMatrix[2][0] += wj * ( lambdas * gradPsi3 * gradPsi1 + sigmaVal * psi3 * psi1 + beta( uh( quadrature[i][0] ) ) * psi3 * psi1 );
      localMatrix[2][1] += wj * ( lambdas * gradPsi3 * gradPsi2 + sigmaVal * psi3 * psi2 + beta( uh( quadrature[i][0] ) ) * psi3 * psi2 );
      localMatrix[2][2] += wj * ( lambdas * gradPsi3 * gradPsi3 + sigmaVal * psi3 * psi3 + beta( uh( quadrature[i][0] ) ) * psi3 * psi3 );
      if ( useTimeDependent )
      {
         double sigmaTimeVal = functions->sigmaTime( currentTime ) / deltaT;
         localMatrix[0][0] += wj * sigmaTimeVal * psi1 * psi1;
         localMatrix[0][1] += wj * sigmaTimeVal * psi1 * psi2;
         localMatrix[0][2] += wj * sigmaTimeVal * psi1 * psi3;
         localMatrix[1][0] += wj * sigmaTimeVal * psi2 * psi1;
         localMatrix[1][1] += wj * sigmaTimeVal * psi2 * psi2;
         localMatrix[1][2] += wj * sigmaTimeVal * psi2 * psi3;
         localMatrix[2][0] += wj * sigmaTimeVal * psi3 * psi1;
         localMatrix[2][1] += wj * sigmaTimeVal * psi3 * psi2;
         localMatrix[2][2] += wj * sigmaTimeVal * psi3 * psi3;
      }
   }
}

void Element::BuildLocalB( )
{
   localB = std::vector<double>( 3, 0.0 );
   double **quadrature = basis.Quadrature( );
   double jacobian = coords[1] - coords[0];
   for ( int i = 0; i < 3; i++ )
   {
      double wj = quadrature[i][1] * jacobian;
      double psi1 = basis.basis1( quadrature[i][0] );
      double psi2 = basis.basis2( quadrature[i][0] );
      double psi3 = basis.basis3( quadrature[i][0] );
      double fVal = f[0] * psi1 + f[1] * psi2 + f[2] * psi3;
      localB[0] += wj * ( fVal * psi1 + theta( uh( quadrature[i][0] ) ) * psi1 + beta( uh( quadrature[i][0] ) ) * uBeta( uh( quadrature[i][0] ) ) * psi1 );
      localB[1] += wj * ( fVal * psi2 + theta( uh( quadrature[i][0] ) ) * psi2 + beta( uh( quadrature[i][0] ) ) * uBeta( uh( quadrature[i][0] ) ) * psi2 );
      localB[2] += wj * ( fVal * psi3 + theta( uh( quadrature[i][0] ) ) * psi3 + beta( uh( quadrature[i][0] ) ) * uBeta( uh( quadrature[i][0] ) ) * psi3 );
      if ( useTimeDependent )
      {
         double sigmaTimeVal = functions->sigmaTime( currentTime ) / deltaT;
         double uhPrev = uhPrevTime( quadrature[i][0] );
         localB[0] += wj * sigmaTimeVal * uhPrev * psi1;
         localB[1] += wj * sigmaTimeVal * uhPrev * psi2;
         localB[2] += wj * sigmaTimeVal * uhPrev * psi3;
      }
   }
}

void Element::BuildLocalMatrixNewton( std::vector<double> q_0 )
{
   localMatrix = std::vector<std::vector<double>>( 3, std::vector<double>( 3, 0.0 ) );

   double** quadrature = basis.Quadrature();
   double jacobian = coords[1] - coords[0];
   for ( int g = 0; g < 3; g++ )
   {
      std::vector<double> psi( 3, 0.0 );
      std::vector<double> gradPsi( 3, 0.0 );
      double wj = quadrature[g][1] * jacobian;

      psi[0] = basis.basis1( quadrature[g][0] );
      psi[1] = basis.basis2( quadrature[g][0] );
      psi[2] = basis.basis3( quadrature[g][0] );

      gradPsi[0] = basis.basis1Grad( quadrature[g][0] ) / jacobian;
      gradPsi[1] = basis.basis2Grad( quadrature[g][0] ) / jacobian;
      gradPsi[2] = basis.basis3Grad( quadrature[g][0] ) / jacobian;

      double uhValue = uh( quadrature[g][0] );
      double sigmaVal = useTimeDependent ? 0.0 : functions->sigma( uhValue );
      double lambdaValue = lambdaApprox( quadrature[g][0] );
      double betaValue = beta( uhValue );
      double sigmaTimeValue = useTimeDependent ? functions->sigmaTime( currentTime ) / deltaT : 0.0;

      double uhGradValue = 0.0;
      for ( int r = 0; r < 3; r++ )
      {
         uhGradValue += q_0[r] * gradPsi[r];
      }

      for ( int i = 0; i < 3; i++ )
      {
         for ( int j = 0; j < 3; j++ )
         {
            double A_ij = lambdaValue * gradPsi[i] * gradPsi[j] +
               sigmaVal * psi[i] * psi[j] +
               betaValue * psi[i] * psi[j] +
               sigmaTimeValue * psi[i] * psi[j];

            double dA_q = materialDLambda( q_0[j] ) * psi[j] * uhGradValue * gradPsi[i];
            localMatrix[i][j] += wj * ( A_ij + dA_q );
         }
      }
   }
}

void Element::BuildLocalBNewton( std::vector<double> q_0 )
{
   localB = std::vector<double>( 3, 0.0 );
   double** quadrature = basis.Quadrature();
   double jacobian = coords[1] - coords[0];
   for ( int g = 0; g < 3; g++ )
   {
      std::vector<double> psi( 3, 0.0 );
      std::vector<double> gradPsi( 3, 0.0 );
      double wj = quadrature[g][1] * jacobian;

      psi[0] = basis.basis1( quadrature[g][0] );
      psi[1] = basis.basis2( quadrature[g][0] );
      psi[2] = basis.basis3( quadrature[g][0] );

      gradPsi[0] = basis.basis1Grad( quadrature[g][0] ) / jacobian;
      gradPsi[1] = basis.basis2Grad( quadrature[g][0] ) / jacobian;
      gradPsi[2] = basis.basis3Grad( quadrature[g][0] ) / jacobian;

      double uhValue = uh( quadrature[g][0] );
      double thetaValue = theta( uhValue );
      double betaValue = beta( uhValue );
      double uBetaValue = uBeta( uhValue );
      double sigmaTimeValue = useTimeDependent ? functions->sigmaTime( currentTime ) / deltaT : 0.0;
      double uhPrev = useTimeDependent ? uhPrevTime( quadrature[g][0] ) : 0.0;
      double fValue = f[0] * psi[0] + f[1] * psi[1] + f[2] * psi[2];

      double uhGradValue = 0.0;
      double dLambdaValue = 0.0;
      for ( int r = 0; r < 3; r++ )
      {
         uhGradValue += q_0[r] * gradPsi[r];
         dLambdaValue += materialDLambda( q_0[r] ) * q_0[r] * psi[r];
      }

      for ( int i = 0; i < 3; i++ )
      {
         double B_i = fValue * psi[i] +
            thetaValue * psi[i] +
            betaValue * uBetaValue * psi[i] +
            sigmaTimeValue * uhPrev * psi[i];

         double dA_q = dLambdaValue * uhGradValue * gradPsi[i];
         localB[i] += wj * ( B_i + dA_q );
      }
   }
}
