#include "Newton.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

double Newton::Condition(int number, double x)
{
   return functions.boundaryValue(number, x);
}

void Newton::Condition1()
{
   if (elements[0].cond == 1)
   {
      matrix.di[0] = 1;
      b[0] = Condition(1, mesh.meshX[0]);;
      for (int i = 0; i < 2; i++)
      {
         matrix.ggu[1][1] = 0;
         matrix.ggu[2][0] = 0;
      }
   }
   if (elements[elemCount - 1].cond == 1)
   {
      matrix.di[n - 1] = 1;
      b[n - 1] = Condition(1, mesh.meshX[elemCount]);;
      for (int i = 0; i < 2; i++)
      {
         matrix.ggl[n - 1][0] = 0;
         matrix.ggl[n - 1][1] = 0;
      }
   }
}

void Newton::Input()
{
   mesh.InputMeshes();
   std::ifstream file("params.txt");
   if (!file)
   {
      std::cerr << "File open error" << std::endl;
      return;
   }
   file >> eps >> maxIter >> boundaryConditions[0] >> boundaryConditions[1];

   elemCount = mesh.n_x - 1;
   n = 3 * elemCount - elemCount + 1;
   matrix.n = n;
   matrix.ggl.resize(n);
   matrix.ggu.resize(n);
   for (int i = 0; i < n; i++)
   {
      matrix.ggl[i] = std::vector<double>(2, 0.0);
      matrix.ggu[i] = std::vector<double>(2, 0.0);
   }
   matrix.di = std::vector<double>(n, 0.0);

   b.resize(n);
   lambdaNodes.resize(n);
   b0 = b;


   std::vector<double> coords;
   elements.resize(elemCount);
   q = std::vector<double>(n, 0.0);
   if (useTimeDependent)
   {
      functions.currentTime = mesh.meshT[0];
      for (int i = 0; i <= elemCount; i++)
      {
         q[i * 2] = functions.u(mesh.meshX[i]);
         if (i < elemCount)
            q[i * 2 + 1] = functions.u((mesh.meshX[i] + mesh.meshX[i + 1]) / 2.0);
      }
   }
   q[0] = Condition(1, mesh.meshX[0]);
   q[n - 1] = Condition(1, mesh.meshX[elemCount]);
   qPrevTime = q;
   std::vector<double> localB(3, 0.0);
   std::vector<double> localLambdaNodes(3, 0.0);
   for (int i = 0; i < elemCount; i++)
   {
      coords = { mesh.meshX[i], mesh.meshX[i + 1] };
      std::vector<double> localQ(3, 0.0);
      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
         {
            localQ[j] = q[3 * i + j];
            localB[j] = b[3 * i + j];
            localLambdaNodes[j] = lambdaNodes[3 * i + j];
         }
         else
         {
            localQ[j] = q[2 * i + j];
            localB[j] = b[2 * i + j];
            localLambdaNodes[j] = lambdaNodes[2 * i + j];
         }
      }
      elements[i] = Element(coords, 0, localQ, localB, localLambdaNodes, &functions);
   }
   elements[0].cond = boundaryConditions[0];
   elements[elemCount - 1].cond = boundaryConditions[1];
   UpdateLayerData();


}

void Newton::BuildMatrix(std::vector<double> qOld)
{
   matrix.di = std::vector<double>(n, 0.0);
   for (int i = 0; i < n; i++)
   {
      matrix.ggl[i] = std::vector<double>(2, 0.0);
      matrix.ggu[i] = std::vector<double>(2, 0.0);
   }

   for (int i = 0; i < elemCount; i++)
   {
      std::vector<double> localQOld(3, 0.0);
      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
            localQOld[j] = qOld[3 * i + j];
         else
            localQOld[j] = qOld[2 * i + j];
      }
      elements[i].useTimeDependent = useTimeDependent;
      elements[i].deltaT = deltaT;
      elements[i].currentTime = mesh.meshT[timeLayer];
      elements[i].BuildLocalMatrixNewton(localQOld);
      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
            matrix.di[3 * i + j] += elements[i].localMatrix[j][j];
         else
            matrix.di[2 * i + j] += elements[i].localMatrix[j][j];
      }
      matrix.ggl[2 * i + 1][1] += elements[i].localMatrix[1][0];
      matrix.ggl[2 * i + 2][0] += elements[i].localMatrix[2][0];
      matrix.ggl[2 * i + 2][1] += elements[i].localMatrix[2][1];

      matrix.ggu[2 * i + 1][1] += elements[i].localMatrix[0][1];
      matrix.ggu[2 * i + 2][0] += elements[i].localMatrix[0][2];
      matrix.ggu[2 * i + 2][1] += elements[i].localMatrix[1][2];
   }
}

void Newton::BuildB(std::vector<double> qOld)
{
   b = std::vector<double>(n, 0.0);
   for (int i = 0; i < elemCount; i++)
   {
      std::vector<double> localQOld(3, 0.0);
      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
            localQOld[j] = qOld[3 * i + j];
         else
            localQOld[j] = qOld[2 * i + j];
      }
      elements[i].useTimeDependent = useTimeDependent;
      elements[i].deltaT = deltaT;
      elements[i].currentTime = mesh.meshT[timeLayer];
      elements[i].BuildLocalBNewton(localQOld);
      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
            b[3 * i + j] += elements[i].localB[j];
         else
            b[2 * i + j] += elements[i].localB[j];
      }
   }
}

void Newton::SolveIter()
{
   qOLd = q;
   UpdateLayerData();
   BuildMatrix(qOLd);
   BuildB(qOLd);
   Condition1();
   lu = LU(matrix, b);
   lu.calcLU();
   lu.calcY();
   lu.calcQ();
   q = ops.AddVec(ops.MultVecScal(lu.q, w), ops.MultVecScal(qOLd, 1 - w));
   std::vector<double> localQ(3, 0.0);
   for (int i = 0; i < elemCount; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
         {
            localQ[j] = q[3 * i + j];
         }
         else
         {
            localQ[j] = q[2 * i + j];
         }
      }
      elements[i].q = localQ;
   }

}

void Newton::Solve()
{
   if (!useTimeDependent)
   {
      int k = 0;
      SolveIter();
      double residual = CalcResidual();

      while (residual > eps && k < maxIter)
      {
         k++;
         residual = CalcResidual();
         std::cout << k << "\tresidual: " << residual << std::endl;
         for (int i = 0; i < elemCount; i++)
            for (int j = 0; j < 3; j++)
            {
               if (i == 0)
                  elements[i].q[j] = q[3 * i + j];
               else
                  elements[i].q[j] = q[2 * i + j];
            }
         SolveIter();
      }

      return;
   }

   timeLayer = 0;
   functions.currentTime = mesh.meshT[0];
   PrintTimeLayerInfo(0, 0);

   for (int s = 1; s < mesh.n_t; s++)
   {
      timeLayer = s;
      deltaT = mesh.meshT[s] - mesh.meshT[s - 1];
      qPrevTime = q;
      q = qPrevTime;

      for (int i = 0; i < elemCount; i++)
      {
         std::vector<double> localQ(3, 0.0);
         std::vector<double> localPrevQ(3, 0.0);
         for (int j = 0; j < 3; j++)
         {
            if (i == 0)
            {
               localQ[j] = q[3 * i + j];
               localPrevQ[j] = qPrevTime[3 * i + j];
            }
            else
            {
               localQ[j] = q[2 * i + j];
               localPrevQ[j] = qPrevTime[2 * i + j];
            }
         }
         elements[i].q = localQ;
         elements[i].qPrevTime = localPrevQ;
      }

      int k = 0;
      SolveIter();
      double residual = CalcResidual();

      while (residual > eps && k < maxIter)
      {
         k++;
         residual = CalcResidual();
         SolveIter();
      }

      functions.currentTime = mesh.meshT[s];
      PrintTimeLayerInfo(s, k);
   }
}

void Newton::UpdateLayerData()
{
   functions.currentTime = mesh.meshT[timeLayer];
   b.assign(n, 0.0);
   lambdaNodes.assign(n, 0.0);

   for (int i = 0; i <= elemCount; i++)
   {
      b[i * 2] = functions.f(mesh.meshX[i]);
      lambdaNodes[i * 2] = functions.lambda(q[i * 2]);
      if (i < elemCount)
      {
         double middle = (mesh.meshX[i] + mesh.meshX[i + 1]) / 2.0;
         b[i * 2 + 1] = functions.f(middle);
         lambdaNodes[i * 2 + 1] = functions.lambda(q[i * 2 + 1]);
      }
   }

   for (int i = 0; i < elemCount; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
         {
            elements[i].f[j] = b[3 * i + j];
            elements[i].lambdaNodes[j] = lambdaNodes[3 * i + j];
         }
         else
         {
            elements[i].f[j] = b[2 * i + j];
            elements[i].lambdaNodes[j] = lambdaNodes[2 * i + j];
         }
      }
   }
}

double Newton::CalcResidual()
{
   UpdateLayerData();

   Matrix A_nonlin;
   A_nonlin.n = n;
   A_nonlin.di.assign(n, 0.0);
   A_nonlin.ggl.assign(n, std::vector<double>(2, 0.0));
   A_nonlin.ggu.assign(n, std::vector<double>(2, 0.0));

   std::vector<double> b_nonlin(n, 0.0);

   for (int i = 0; i < elemCount; i++)
   {
      elements[i].useTimeDependent = useTimeDependent;
      elements[i].deltaT = deltaT;
      elements[i].currentTime = mesh.meshT[timeLayer];
      elements[i].BuildLocalMatrix();
      elements[i].BuildLocalB();

      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
            A_nonlin.di[3 * i + j] += elements[i].localMatrix[j][j];
         else
            A_nonlin.di[2 * i + j] += elements[i].localMatrix[j][j];
      }

      A_nonlin.ggl[2 * i + 1][1] += elements[i].localMatrix[1][0];
      A_nonlin.ggl[2 * i + 2][0] += elements[i].localMatrix[2][0];
      A_nonlin.ggl[2 * i + 2][1] += elements[i].localMatrix[2][1];

      A_nonlin.ggu[2 * i + 1][1] += elements[i].localMatrix[0][1];
      A_nonlin.ggu[2 * i + 2][0] += elements[i].localMatrix[0][2];
      A_nonlin.ggu[2 * i + 2][1] += elements[i].localMatrix[1][2];

      for (int j = 0; j < 3; j++)
      {
         if (i == 0)
            b_nonlin[3 * i + j] += elements[i].localB[j];
         else
            b_nonlin[2 * i + j] += elements[i].localB[j];
      }
   }

   if (elements[0].cond == 1)
   {
      A_nonlin.di[0] = 1;
      b_nonlin[0] = Condition(1, mesh.meshX[0]);
      A_nonlin.ggu[1][1] = 0;
      A_nonlin.ggu[2][0] = 0;
   }
   if (elements[elemCount - 1].cond == 1)
   {
      A_nonlin.di[n - 1] = 1;
      b_nonlin[n - 1] = Condition(1, mesh.meshX[elemCount]);
      A_nonlin.ggl[n - 1][0] = 0;
      A_nonlin.ggl[n - 1][1] = 0;
   }

   std::vector<double> r = ops.AddVec(
      ops.MultMatVec(A_nonlin, q),
      ops.MultVecScal(b_nonlin, -1)
   );

   return ops.DotProduct(r, r) / ops.DotProduct(b_nonlin, b_nonlin);
}

double Newton::CalcL2Error()
{
   double sum = 0.0;
   for (int i = 0; i < elemCount; i++)
   {
      double** quadrature = elements[i].basis.Quadrature();
      double jacobian = elements[i].coords[1] - elements[i].coords[0];
      for (int g = 0; g < 3; g++)
      {
         double xi = quadrature[g][0];
         double weight = quadrature[g][1] * jacobian;
         double x = elements[i].coords[0] + xi * jacobian;
         double diff = elements[i].uh(xi) - functions.u(x);
         sum += weight * diff * diff;
      }
   }
   return std::sqrt(sum);
}

void Newton::PrintTimeLayerInfo(int layerIndex, int iterationCount)
{
   double error = CalcL2Error();
   std::streamsize oldPrecision = std::cout.precision();
   std::ios::fmtflags oldFlags = std::cout.flags();

   std::cout << std::fixed << std::setprecision(5)
      << "Time layer " << layerIndex << ": t = " << mesh.meshT[layerIndex]
      << ", countNodes = " << mesh.n_x;

   if (layerIndex > 0)
   {
      std::cout << ", ht = " << (mesh.meshT[layerIndex] - mesh.meshT[layerIndex - 1])
         << ", countIteration = " << iterationCount;
   }

   std::cout << std::scientific << std::setprecision(4)
      << ", ||u* - u||L2 = " << error << std::endl;

   std::cout.flags(oldFlags);
   std::cout.precision(oldPrecision);
}

