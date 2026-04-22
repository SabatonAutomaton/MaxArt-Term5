#include "PenaltyBarrierMethod.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

static fs::path gProjectRoot;

static fs::path FindProjectRoot()
{
   fs::path cur = fs::current_path();
   while (!cur.empty())
   {
      if (fs::exists(cur / "in.txt") && fs::exists(cur / "MO_3.vcxproj"))
         return cur;
      fs::path parent = cur.parent_path();
      if (parent == cur)
         break;
      cur = parent;
   }
   return fs::current_path();
}

static const char* MethodName(int methodType)
{
   return (methodType == 1) ? "penalty" : "barrier";
}

static const char* ProblemName(int methodType, int problemType)
{
   if (methodType == 1)
      return (problemType == 1) ? "task_a" : "task_b";
   return "task_a";
}

static const char* ModifierName(int methodType, int variant)
{
   if (methodType == 1)
   {
      switch (variant)
      {
      case 1: return "max(0,g)";
      case 2: return "max(0,g)^2";
      case 3: return "max(0,g)^4";
      default: return "max(0,g)^2";
      }
   }

   switch (variant)
   {
   case 1: return "-ln(-g)";
   case 2: return "-1/g";
   case 3: return "1/g^2";
   default: return "-ln(-g)";
   }
}

static const char* StrategyName(int methodType, int variant)
{
   if (methodType == 1)
   {
      switch (variant)
      {
      case 1: return "r_{k+1}=beta*r_k";
      case 2: return "r_{k+1}=r_k+beta";
      case 3: return "r_{k+1}=r_k+1";
      default: return "r_{k+1}=beta*r_k";
      }
   }

   switch (variant)
   {
   case 1: return "r_{k+1}=r_k/beta";
   case 2: return "r_{k+1}=max(r_k-beta,eps)";
   case 3: return "r_{k+1}=r_k/2";
   default: return "r_{k+1}=r_k/beta";
   }
}

static void WriteSummaryHeaderIfNeeded()
{
   const fs::path summaryFile = gProjectRoot / "out1.txt";
   bool rewrite = true;
   if (fs::exists(summaryFile))
   {
      std::ifstream in(summaryFile);
      std::string firstLine;
      std::getline(in, firstLine);
      rewrite = firstLine.find("strategy_r") == std::string::npos || firstLine.find("func_evals") == std::string::npos;
   }

   if (rewrite)
   {
      std::ofstream out(summaryFile, std::ios::trunc);
      out << "method\tproblem\tmodifier\tstrategy_r\tstart_x\tstart_y\teps\tr0\tbeta\titers\tinner_iters\tfunc_evals\tresult_x\tresult_y\tobjective\tviolation_or_gap\tcoefficient" << std::endl;
   }
}

static void AppendSummary(int methodType, int problemType, int modifierVariant, int strategyVariant, const double* startPoint,
   double eps, double r0, double beta, int outerIters, int innerIters, long long funcEvals,
   double resultX, double resultY, double objective, double measure, double coefficient)
{
   WriteSummaryHeaderIfNeeded();
   std::ofstream out(gProjectRoot / "out1.txt", std::ios::app);
   out << MethodName(methodType) << '\t'
      << ProblemName(methodType, problemType) << '\t'
      << ModifierName(methodType, modifierVariant) << '\t'
      << StrategyName(methodType, strategyVariant) << '\t'
      << std::setprecision(16) << startPoint[0] << '\t'
      << startPoint[1] << '\t'
      << eps << '\t'
      << r0 << '\t'
      << beta << '\t'
      << outerIters << '\t'
      << innerIters << '\t'
      << funcEvals << '\t'
      << resultX << '\t'
      << resultY << '\t'
      << objective << '\t'
      << measure << '\t'
      << coefficient << std::endl;
}

static void RunPlotter()
{
   fs::path script = gProjectRoot.parent_path() / "Graphiks_MO3" / "Graphiks_MO3.py";
   if (!fs::exists(script))
      return;

   std::string command = "py \"" + script.string() + "\"";
   if (std::system(command.c_str()) != 0)
   {
      command = "python \"" + script.string() + "\"";
      std::system(command.c_str());
   }
}

int main()
{
   gProjectRoot = FindProjectRoot();
   fs::current_path(gProjectRoot);

   PenaltyBarrierMethod solver;
   int methodType;
   int problemType = 1;
   int modifierVariant = 2;
   int strategyVariant = 1;
   double startPoint[2];
   double eps;
   double r0 = 1.0;
   double beta = 10.0;
   std::ifstream input(gProjectRoot / "in.txt");

   if (!input)
   {
      std::cerr << "Cannot open in.txt" << std::endl;
      return 1;
   }

   if (!(input >> startPoint[0] >> startPoint[1] >> eps))
   {
      std::cerr << "Invalid format in.txt. Expected: x0 y0 eps [r0 beta]" << std::endl;
      return 1;
   }

   if (!(input >> r0 >> beta))
   {
      r0 = 1.0;
      beta = 10.0;
   }

   std::cout << "Method (1 - penalty, 2 - barrier): ";
   std::cin >> methodType;

   if (methodType == 1)
   {
      std::cout << "Penalty problem (1 - task a: y - x >= 1, 2 - task b: x = -y): ";
      std::cin >> problemType;
      std::cout << "Penalty function (1 - max(0,g), 2 - max(0,g)^2, 3 - max(0,g)^4): ";
      std::cin >> modifierVariant;
      std::cout << "Penalty r strategy (1 - beta*r_k, 2 - r_k+beta, 3 - r_k+1): ";
      std::cin >> strategyVariant;

      solver.InitPenalty(2, problemType);
      solver.SetModifierVariant(modifierVariant);
      solver.SetStrategyVariant(strategyVariant);
      solver.SetStartPoint(startPoint);
      solver.SetTolerance(eps);
      solver.SetParameters(r0, beta);
   }
   else
   {
      problemType = 1;
      std::cout << "Barrier function (1 - -ln(-g), 2 - -1/g, 3 - 1/g^2): ";
      std::cin >> modifierVariant;
      std::cout << "Barrier r strategy (1 - r_k/beta, 2 - max(r_k-beta,eps), 3 - r_k/2): ";
      std::cin >> strategyVariant;

      solver.InitBarrier(2);
      solver.SetModifierVariant(modifierVariant);
      solver.SetStrategyVariant(strategyVariant);
      solver.SetStartPoint(startPoint);
      solver.SetTolerance(eps);
      solver.SetParameters(r0, beta);
   }

   double* res = solver.Solver();

   std::cout << "\nResult: ";
   for (int i = 0; i < 2; ++i)
      std::cout << std::setprecision(16) << res[i] << " ";
   std::cout << "\nOuter iterations: " << solver.GetLastOuterIterationCount();
   std::cout << "\nTotal inner iterations: " << solver.GetTotalInnerIterationCount();
   std::cout << "\nObjective value: " << std::setprecision(16) << solver.GetLastObjectiveValue();
   std::cout << "\n" << ((methodType == 1) ? "Constraint violation: " : "Boundary gap: ")
      << std::setprecision(16) << solver.GetLastMeasure();
   std::cout << "\nLast coefficient: " << std::setprecision(16) << solver.GetLastCoefficient();
   std::cout << "\nFunction evaluations: " << solver.GetTotalFunctionEvaluationCount() << std::endl;

   AppendSummary(methodType, problemType, modifierVariant, strategyVariant, startPoint, eps, r0, beta,
      solver.GetLastOuterIterationCount(), solver.GetTotalInnerIterationCount(), solver.GetTotalFunctionEvaluationCount(),
      res[0], res[1], solver.GetLastObjectiveValue(), solver.GetLastMeasure(), solver.GetLastCoefficient());

   RunPlotter();
   return 0;
}
