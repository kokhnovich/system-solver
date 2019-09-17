#include <bits/stdc++.h>
using namespace std;

template<typename T>
using Matrix = vector<vector<T>>;

template<typename T>
vector<vector<T>> mult(vector<vector<T>> A, vector<vector<T>> B) {
  if (A.empty() || B.empty() || A[0].size() != B.size()) throw logic_error("wrong sizes");
  size_t N = A.size(), M = A[0].size(), K = B[0].size();
  vector<vector<T>> R(N, vector<T>(K, T(0)));
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < K; j++) {
      for (size_t k = 0; k < M; k++) {
        R[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return R;
}

// @TODO add extended matrix solver method
// @TODO norm tests
// @TODO optimize using profiler
template<class T>
struct greater_using_abs {
  bool operator()(const T& x, const T& y) const { return abs(x) > abs(y); }
};

template<typename T>
vector<vector<T>> getHWMatrix(int n) {
  return {{n + 1, n / 2, -n / 2, 1},
          {-n - 1, -(n + 1) / 2, (n + 1) / 2, -(n + 2) / 3},
          {-n + 1, (n + 1) / 2, -(n + 2) / 3, n - 1},
          {n / 3, -1, n, -n}};
}

template<typename T>
Matrix<T> getIdentityMatrix(int n) {
  vector<vector<T>> A(n, vector<T>(n, T(0)));
  for (size_t i = 0; i < n; ++i) {
    A[i][i] = T(1);
  }
  return A;
}

#include "Fraction.h"
#include "Solver.h"
#include "Debug.h"
#include "test_runner.h"

void SolveMyHomework(int n = 12) {
  /*
  cout << "n = 12\nmatrix is\n";
  vector<vector<Fraction>> a = getHWMatrix<Fraction>(n);
  vector<Fraction> b(4);
  for (size_t i = 0; i < b.size(); ++i) {
    // b[i] = accumulate(a[i].begin(), a[i].end(), Fraction(0, 1));
    b[i] = a[i][0] + Fraction(1) * a[i][1] + Fraction(1) * a[i][2] + Fraction(1) * a[i][3];
  }
  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>();
  for (auto& i : solver->SolveSystem(a, b, SolverMethod::BEST_IN_MATRIX)) {
    cout << i << " ";
  }
  cout << endl;
   */
}

template<typename T>
T cubeNorm(const vector<vector<T>>& A) {
  T best = T(0);
  for (int j = 0; j < A[0].size(); ++j) {
    T now = T(0);
    for (int i = 0; i < A.size(); ++i) {
      now += abs(A[i][j]);
    }
    best = max(best, now);
  }
  return best;
}
/*
void SolveMyHomeWorkAboutLUP() {

  auto A = getHWMatrix<Fraction>(12);
  vector<vector<Fraction>> L, U, P;
  vector<Fraction> b(4);
  for (size_t i = 0; i < b.size(); ++i) {
    // b[i] = accumulate(a[i].begin(), a[i].end(), Fraction(0, 1));
    b[i] = A[i][0] + Fraction(1) * A[i][1] + Fraction(1) * A[i][2] + Fraction(1) * A[i][3];
  }

  PrintMatrix(A);
  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>();
  tie(L, U, P) = solver->LUP_Decomposition(A, SolverMethod::BEST_IN_ROW);

  PrintMatrix(L, "L");

  PrintMatrix(U, "U");

  PrintMatrix(P, "P");

  auto res = mult(L, U);
  PrintMatrix(res, "LU");
  PrintMatrix(mult(U, P), "UP");
  res = mult(res, P);
  PrintMatrix(res, "A = LUP");
  PrintMatrix(mult(L, mult(U, P)), "A = L(UP)");

  cout << endl << (A == res ? "works nice (:" : "smth goes wrong ):") << endl;

  auto B = getIdentityMatrix<Fraction>(4);
  auto ans = solver->SolveSystemUsingLU(A, B, SolverMethod::BEST_IN_ROW);

  PrintMatrix(ans, "X");
  PrintMatrix(mult(A, ans), "AX");

  cout << "Cube norm " << cubeNorm(ans) * Fraction(41, 1) << endl;

  cout << endl << (B == mult(A, ans) ? "good job bro" : "stupid wood") << endl;

  ans = solver->SolveSystem(A, B, SolverMethod::BEST_IN_MATRIX);

  PrintMatrix(ans, "X 2.0");
  PrintMatrix(mult(A, ans), "AX 2.0");

  cout << endl << (B == mult(A, ans) ? "good job bro" : "stupid wood") << endl;
}
*/

void SolveMyZALUPHomeWorkAboutDLUP() {

  Matrix<Fraction> A = {{1, 3, 2}, {4, 5, 4}, {3, 5, 3}}, L, D, U, P;

  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>();

  //PrintMatrix(solver->GetReversed(A), "A^{-1}");

  PrintMatrix(A, "A before");
  tie(D, L, U, P) = solver->DLUP_Decomposition(A);
  PrintMatrix(A, "A after");
  PrintMatrix(D, "D");
  PrintMatrix(L, "L");
  PrintMatrix(U, "U");
  PrintMatrix(P, "P");

  PrintMatrix(mult(mult(L, D), mult(U, P)), "LDUP");
  PrintMatrix(mult(D, mult(U, P)), "DUP");
  PrintMatrix(mult(L, U), "LU");
  PrintMatrix(mult(mult(D, A), P), "DAP");
}

int main() {
//  // SolveMyHomework();
//  vector<int> a(2, 3), b(2, 2);
//
//  // cout << mult<int>({{1, 1}, {1, 1}, {1, 1}}, {{1, 1, 1}, {1, 1, 1}}) << endl;
//
//  cout << getIdentityMatrix<Fraction>(4) << endl;
//
//  vector<vector<Fraction>> A = {{1, 2}, {3, 4}}, B{{1, 0}, {0, 1}};
//  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>();
//  auto ans = solver->SolveSystem(A, B, SolverMethod::DO_NOT_TOUCH);
//  cout << ans.size() << " " << ans[1].size() << endl;
//  cout << mult(A, ans) << endl;
//  for (const auto& i : ans) {
//    for (const auto& j : i) {
//      cout << j << ' ';
//    }
//    cout << endl;
//  }

  SolveMyZALUPHomeWorkAboutDLUP();
  return 0;
}