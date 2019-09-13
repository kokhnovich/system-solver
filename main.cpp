#include <bits/stdc++.h>
#include "profile.h"
#include "test_runner.h"
#include "Fraction.h"

using namespace std;

enum class SolverMethod {
  DO_NOT_TOUCH,
  BEST_IN_ROW,
  BEST_IN_COLUMN,
  BEST_IN_MATRIX
};

template<typename T>
vector<vector<T>> mult(vector<vector<T>> A, vector<vector<T>> B) {
  size_t n = A.size();
  vector<vector<T>> R(n, vector<T>(n, T(0)));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < n; k++) {
        R[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return R;
}

/** matrix should be NxN
 *  template class Func can be, for example, std:greater<T> or std:less<T>
 *  Func must be linear!
 *
**/

// @TODO file separation
// @TODO norm tests
// @TODO optimize using profiler
// @TODO add LU-decomposition
// @TODO architecture + memory optimizations


Fraction abs(const Fraction& val) {
  return Fraction(abs(val.numerator), abs(val.denominator));
}

template<class T>
struct greater_using_abs {
  bool operator()(const T& x, const T& y) const { return abs(x) > abs(y); }
};

#include "Solver.h"

template<typename T>
void PrintMatrix(const vector<vector<T>>& a, const string& message = "") {
  cout << message << ":\n";
  for (auto& i : a) {
    for (auto& j : i) {
      cout << j << " ";
    }
    cout << endl;
  }
  cout << "\n\n";
}

template<typename T>
vector<vector<T>> getHWMatrix(int n) {
  return {{n + 1, n / 2, -n / 2, 1},
          {-n - 1, -(n + 1) / 2, (n + 1) / 2, -(n + 2) / 3},
          {-n + 1, (n + 1) / 2, -(n + 2) / 3, n - 1},
          {n / 3, -1, n, -n}};
}

void SolveMyHomework(int n = 12) {
  cout << "n = 12\nmatrix is\n";
  vector<vector<Fraction>> a = getHWMatrix<Fraction>(n);
  vector<Fraction> b(4);
  for (size_t i = 0; i < b.size(); ++i) {
    // b[i] = accumulate(a[i].begin(), a[i].end(), Fraction(0, 1));
    b[i] = a[i][0] + Fraction(1) * a[i][1] + Fraction(1) * a[i][2] + Fraction(1) * a[i][3];
  }
  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>(a, b, SolverMethod::BEST_IN_MATRIX);
  solver->Print();
  for (auto& i : solver->GetSolution()) {
    cout << i << " ";
  }
  cout << endl;
  solver->Print();
}

void SolveMyHomeWorkAboutLUP() {
  auto A = getHWMatrix<Fraction>(12);
  vector<vector<Fraction>> L, U, P;
  vector<Fraction> b(4);
  for (size_t i = 0; i < b.size(); ++i) {
    // b[i] = accumulate(a[i].begin(), a[i].end(), Fraction(0, 1));
    b[i] = A[i][0] + Fraction(1) * A[i][1] + Fraction(1) * A[i][2] + Fraction(1) * A[i][3];
  }

  PrintMatrix(A);
  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>(A, b, SolverMethod::DO_NOT_TOUCH);
  tie(L, U, P) = solver->LUP_Decomposition();

  PrintMatrix(L, "L");

  PrintMatrix(U, "U");

  PrintMatrix(P, "P");

  auto res = mult(L, U);
  PrintMatrix(res, "LU");
  PrintMatrix(mult(U, P), "UP");
  res = mult(res, P);
  PrintMatrix(res, "A = LUP");

  cout << endl << (A == res ? "works nice (:" : "smth goes wrong ):") << endl;

  auto ans = solver->SolveSystemUsingLU();

  cout << "ans\n";
  for (const auto& i : ans) {
    cout << i << " ";
  }
  cout << endl;

}

#include "tests.cpp"

int main() {
  // SolveMyHomework();
  SolveMyHomeWorkAboutLUP();
  return 0;
}