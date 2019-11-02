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

// #define double long double

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
// @TODO iteration methods
// @TODO make directories

#include "utils.cpp"
#include "Fraction.h"
#include "test_runner.h"
#include "Debug.h"
#include "Solver.cpp"
#include "tests.cpp"
#include "profile.h"
#include "hwsolver.h"
#include "LatexWriter.h"
#include "IterationMethods/IterationMethods.h"
#include "lab/cma_lab_tasks.cpp"

#include "MatrixChanger.h"

int main() {
  // SolveTask1();
  SolveTask2();
  //SolveTask3();
  // SolveTask4();
  // SolveTask5();
  return 0;
}
