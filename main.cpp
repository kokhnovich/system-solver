//#include <iostream>
//#include <vector>
//#include <random>
//#include <cassert>
//#include "algorithm"
//#include <iomanip>
//#include <chrono>
//#include <future>
//#include <functional>
#include <bits/stdc++.h>

// #pragma comment(linker, "/stack:200000000")
// #pragma GCC optimize("O3")
// #pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

using namespace std;
using namespace placeholders;

template<typename T>
using Matrix = vector<vector<T>>;


// @TODO boost storage of rows/column-matrixes
// @TODO LABAAA
// @TODO norm tests
// @TODO optimize using profiler
// @TODO class LatexWriter


#include "utils.cpp"
#include "Fraction.h"
#include "test_runner.h"
#include "Debug.h"
#include "Solver.cpp"
#include "tests.cpp"
#include "profile.h"
#include "hwsolver.h"
#include "cma_lab_tasks.cpp"

int main() {
//  PrintMatrix(generateRandomSymmetricMatrix<Fraction>(10), "random sym matrix");
//  // SolveMyHomeWorkAboutDLUP();
//  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>();
//  vector<ThreeDiagonal<Fraction>>
//      a = {{0, 2, 3},
//           {1, 3, 2},
//           {2, 2, 8},
//           {4, 2, 4},
//           {2, 4, 0}};
//  vector<Fraction> b = {5, 6, 12, 10, 6};
//  solver->SolveThreeDiagonalSystem(a, b);

  // SolveTask1();
  // SolveTask2();
  // SolveTask3();
  SolveTask4();
  return 0;
}
