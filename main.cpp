#include <bits/stdc++.h>
using namespace std;

template<typename T>
using Matrix = vector<vector<T>>;

template<typename T>
Matrix<T> mult(Matrix<T> A, Matrix<T> B) {
  if (A.empty() || B.empty() || A[0].size() != B.size()) throw logic_error("wrong sizes");
  size_t N = A.size(), M = A[0].size(), K = B[0].size();
  Matrix<T> R(N, vector<T>(K, T(0)));
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < K; j++) {
      for (size_t k = 0; k < M; k++) {
        R[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return R;
}

// @TODO norm tests
// @TODO optimize using profiler
// @TODO class LatexWriter

template<class T>
struct greater_using_abs {
  bool operator()(const T& x, const T& y) const { return abs(x) > abs(y); }
};

template<typename T>
Matrix<T> getHWMatrix(int n) {
  return {{n + 1, n / 2, -n / 2, 1},
          {-n - 1, -(n + 1) / 2, (n + 1) / 2, -(n + 2) / 3},
          {-n + 1, (n + 1) / 2, -(n + 2) / 3, n - 1},
          {n / 3, -1, n, -n}};
}

template<typename T>
Matrix<T> getIdentityMatrix(int n) {
  Matrix<T> A(n, vector<T>(n, T(0)));
  for (size_t i = 0; i < n; ++i) {
    A[i][i] = T(1);
  }
  return A;
}

#include "Fraction.h"
#include "test_runner.h"
#include "Solver.cpp"
#include "Debug.h"
#include "tests.cpp"

void SolveMyHomeWorkAboutDLUP() {
  Matrix<Fraction> A = {{1, 3, 2},
                        {3, 5, 7},
                        {4, 5, 8}}, L, D, U, P;
  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>();

  PrintMatrix(A, "A before");
  tie(D, L, U, P) = solver->DLUP_Decomposition(A);
  PrintMatrix(A, "A after");
  PrintMatrix(D, "D");
  PrintMatrix(L, "L");
  PrintMatrix(U, "U");
  PrintMatrix(P, "P");

  auto dlup = mult(solver->GetReversed(D), mult(L, mult(U, solver->GetReversed(P))));
  PrintMatrix(dlup, "D^{-1}LUP^{-1}");

  dlup = mult(mult(solver->GetReversed(D), L), mult(U, solver->GetReversed(P)));
  PrintMatrix(dlup, "(D^{-1}L) (UP^{-1})");

  cout << (dlup == A ? "Wonderful!" : "wrong") << endl;
  assert(dlup == A);
  PrintMatrix(mult(L, U), "LU");
  PrintMatrix(mult(mult((D), A), (P)), "DAP");

  PrintMatrix(mult(D, P), "DP");
}

int main() {
  // SolveMyHomeWorkAboutDLUP();
  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>();
  vector<ThreeDiagonal<Fraction>> a = {{0, 2, 3},
                                       {1, 3, 2},
                                       {2, 2, 8},
                                       {4, 2, 4},
                                       {2, 4, 0}};
  vector<Fraction> b = {5, 6, 12, 10, 6};
  solver->SolveThreeDiagonalSystem(a, b);
  return 0;
}