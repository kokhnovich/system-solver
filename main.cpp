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
// @TODO iteration methods

#include "utils.cpp"
#include "Fraction.h"
#include "test_runner.h"
#include "Debug.h"
#include "Solver.cpp"
#include "tests.cpp"
#include "profile.h"
//#include "hwsolver.h"
//#include "LatexWriter.h"
//#include "IterationMethods/IterationMethods.h"
//#include "cma_lab_tasks.cpp"

#include "MatrixChanger.h"

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

  // LatexWriter latex_writer("latex.tex");
  // latex_writer.Testing();

  // SolveTask1();
  // SolveTask2();
  // SolveTask3();
  // SolveTask4();
  // SolveTask5();

//  Matrix<double> test = {{1, 2, 3},
//                         {4, 5, 6},
//                         {7, 8, 9}};
//  auto changer = MatrixChanger(test);
//
//  changer.Print();
//
//  changer.MultRow(1, 25);
//  changer.Print();
//
//  changer.AddRow(1, 0, 1);
//
//  changer.Print();
//
//  return 0;

  auto to_solve = getNewHWMatrix<double>(13);

  Matrix<double> jopa = {{-11, 12},
                         {-24, 25}};

  MatrixChanger matrix_changer(jopa);

  matrix_changer.Print();
  while (true) {
    string cmd;
    cin >> cmd;
    if (cmd == "swaprow") {
      int i, j;
      cin >> i >> j;
      matrix_changer.SwapRows(i, j);
    } else if (cmd == "addrow") {
      int i, j;
      double k;
      cin >> i >> j >> k;
      matrix_changer.AddRow(i, j, double(k));
    } else if (cmd == "multrow") {
      int i;
      double k;
      cin >> i >> k;
      matrix_changer.MultRow(i, double(k));
    } else if (cmd == "swapcol") {
      int i, j;
      cin >> i >> j;
      matrix_changer.SwapCols(i, j);
    } else if (cmd == "addcol") {
      int i, j;
      double k;
      cin >> i >> j >> k;
      matrix_changer.AddCol(i, j, double(k));
    } else if (cmd == "multcol") {
      int i;
      double k;
      cin >> i >> k;
      matrix_changer.MultCol(i, double(k));
    }
    matrix_changer.Print();
  }

  /*
  IterationMethods<double> solver;

  Matrix<double> A = {{200, 4, -4},
                      {-3, 150, 3},
                      {2, -2, 100}};

  Matrix<double> x = {{1}, {-2}, {1}};

  auto b = mult(A, x);
  PrintMatrix(b, "b");
  PrintMatrix(solver.YakobyMethod(A, make1Dfrom2D(b)), "yakoby");
  PrintMatrix(solver.GaussZeydelMethod(A, make1Dfrom2D(b)), "gauss");
  PrintMatrix(solver.RelaxationMethod(A, make1Dfrom2D(b), 0.9), "relax");
  */

  return 0;
}
