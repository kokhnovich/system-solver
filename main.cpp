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

template<typename T>
Matrix<T> generateRandomSymmetricMatrix(int size_, int min_elem = -100, int max_elem = 100) {
  std::uniform_int_distribution<std::mt19937_64::result_type> udist(min_elem, max_elem);
  std::mt19937_64 generator(std::random_device{}());
  Matrix<T> a(size_, vector<T>(size_));
  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < size_; ++j) {
      if (j <= i) {
        a[i][j] = a[j][i];
      } else {
        a[i][j] = udist(generator);
      }
    }
  }
  return a;
}

template<typename T>
vector<T> generateAnsMatrix(const Matrix<T>& A, const vector<T>& x) {
  vector<T> ans(A.size());
  for (int i = 0; i < A.size(); ++i) {
    ans[i] = T(0);
    for (int j = 0; j < A.size(); ++j) {
      ans[i] += x[j] * A[i][j];
    }
  }
  return ans;
}

/// Nx1 == true  <=> ans.size() == n && ans[i].size() == 1 \forall i
/// Nx1 == false <=> ans.size() == 1 && ans[i].size() == n \forall i
template<typename T>
Matrix<T> make2Dfrom1D(const vector<T>& a, bool Nx1 = true) {
  if (Nx1) {
    Matrix<T> ans(a.size(), vector<T>(1));
    for (int i = 0; i < a.size(); ++i) {
      ans[i][0] = a[i];
    }
    return ans;
  } else {
    Matrix<T> ans(1, vector<T>(a.size()));
    for (int i = 0; i < a.size(); ++i) {
      ans[0][i] = a[i];
    }
    return ans;
  }
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
#include "profile.h"

void SolveMyHomeWorkAboutDLUP() {
  Matrix < Fraction > A = {{1, 3, 2},
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

#include "cma_lab_tasks.cpp"
int main() {
  PrintMatrix(generateRandomSymmetricMatrix<Fraction>(10), "random sym matrix");
  // SolveMyHomeWorkAboutDLUP();
  auto solver = new Solver<Fraction, greater_using_abs<Fraction>>();
  vector<ThreeDiagonal<Fraction>>
      a = {{0, 2, 3},
           {1, 3, 2},
           {2, 2, 8},
           {4, 2, 4},
           {2, 4, 0}};
  vector<Fraction> b = {5, 6, 12, 10, 6};
  solver->SolveThreeDiagonalSystem(a, b);
  return 0;
}