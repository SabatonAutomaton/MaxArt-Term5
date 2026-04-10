#include "Problem.h"
#include "LU.h"
#include "MSG.h"
#include "SLE.h"
#include "SLEAssistant.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>


int main( )
{
   int testCase = 1;
   std::ifstream testFile( "test.txt" );
   if ( testFile )
   {
      testFile >> testCase;
   }


   switch ( testCase )
   {
   case 1:
      std::cout << "Running test case 1: 1D manufactured test\n";
      break;
   case 2:
      std::cout << "Running test case 2: 3D harmonic linear manufactured test #1\n";
      break;
   case 3:
      std::cout << "Running test case 3: 3D harmonic linear manufactured test #2\n";
      break;
   default:
      std::cout << "Unknown test in test.txt, fallback to case 1\n";
      testCase = 1;
      break;
   }

   Problem problem;
   problem.Run( testCase );
   return 0;
}
