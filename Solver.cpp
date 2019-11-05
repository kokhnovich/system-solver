//
// Created by user on 25.09.19.
//
#include "Solver.h"

template<typename T, class Func>
tuple<Matrix<T>, Matrix<T>, Matrix<T>> Solver<T, Func>::LUP_Decomposition(const Matrix<T>& A,
                                                                          const SolverMethod& method_) {
  vector<int> ans_order(A.size());
  Matrix<T> L(A.size(), vector<T>(A.size(), 0)), U(A), P(A.size(), vector<T>(A.size(), 0));
  for (int i = 0; i < A.size(); ++i) {
    ans_order[i] = i;
  }
  // Print(U, ans_order);
  for (int stage = 0; stage < U.size(); ++stage) {
    L[stage][stage] = 1;
    switch (method_) {
      case SolverMethod::BEST_IN_ROW : {
        swap_columns(U, stage, best_in_the_row(U, stage), ans_order);
        break;
      }
      case SolverMethod::DO_NOT_TOUCH: {
        break;
      }
      default: {
        throw logic_error("TBD");
      }
    }
    if (U[stage][stage] == T(0)) continue;

    for (int row = stage + 1; row < A.size(); ++row) {
      T k = U[row][stage] / U[stage][stage];
      L[row][stage] = k;
      sub_row(U, row, stage, k);
    }
  }

  for (int i = 0; i < A.size(); ++i) {
    P[i][ans_order[i]] = 1;
  }

  // assert(compareMatrixOfDouble(mult(mult(L, U), P), A));
  return make_tuple(L, U, P);
}

template<typename T, class Func>
void Solver<T, Func>::substract_str(Matrix<T>& A, int row, int stage) {
  T koef = A[row][stage];
  for (int col = stage; col < A[0].size(); ++col) {
    A[row][col] -= A[stage][col] * koef;
  }
}

template<typename T, class Func>
Matrix<T> Solver<T, Func>::SolveSystemUsingLU(Matrix<T> A,
                                              Matrix<T> B,
                                              const SolverMethod& method_) {
  Matrix<T> L, U, P;
  tie(L, U, P) = LUP_Decomposition(A, method_);
//  PrintMatrix(L, "L");
//  PrintMatrix(U, "U");
//  PrintMatrix(P, "P");
  B = LU_Step(B, L, U);
  return B;
}
template<typename T, class Func>
Matrix<T>& Solver<T, Func>::LU_Step(Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& U) {
  for (int i = 0; i < B.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      sub_row(B, i, j, L[i][j]);
    }
    if (L[i][i] != T(1)) {
      for (auto& j : B[i]) {
        j /= L[i][i];
      }
    }
  }
  for (int i = B.size() - 1; i >= 0; --i) {
    for (int j = i + 1; j < B.size(); ++j) {
      sub_row(B, i, j, U[i][j]);
    }
    if (U[i][i] != T(1)) {
      for (auto& j : B[i]) {
        j /= U[i][i];
      }
    }
  }
  return B;
}

template<typename T, class Func>
int Solver<T, Func>::best_in_the_row(const Matrix<T>& A, int row) {
  int best_column = row;
  for (int column = row + 1; column < A[0].size(); ++column) {
    if (Func()(A[row][column], A[row][best_column])) {
      best_column = column;
    }
  }
  return best_column;
}
template<typename T, class Func>
int Solver<T, Func>::best_in_the_col(const Matrix<T>& A, int col) {
  int best_row = col;
  for (int i = col + 1; i < A.size(); ++i) {
    if (Func()(A[i][col], A[best_row][col])) {
      best_row = i;
    }
  }
  return best_row;
}
template<typename T, class Func>
pair<int, int> Solver<T, Func>::best_in_the_sqr(const Matrix<T>& A, int start_i, int start_j) {
  int best_i = start_i, best_j = start_j;
  for (int i = start_i; i < A.size(); ++i) {
    for (int j = start_j; j < A.size(); ++j) {
      if (Func()(A[i][j], A[best_i][best_j])) {
        best_i = i, best_j = j;
      }
    }
  }
  return make_pair(best_i, best_j);
}
template<typename T, class Func>
void Solver<T, Func>::Print(const Matrix<T>& A, const vector<int>& ans_order) {
  for (int i = 0; i < A.size(); ++i) {
    for (auto& j : A[i]) {
      cout << j << " ";
    }
    cout << endl;
  }
  cout << "Порядок неизвестных: ";
  for (const auto& i : ans_order) {
    cout << i + 1 << " ";
  }
  cout << endl;
}
template<typename T, class Func>
void Solver<T, Func>::Print(const Matrix<T>& A) {
  for (int i = 0; i < A.size(); ++i) {
    for (auto& j : A[i]) {
      cout << j << " ";
    }
    cout << endl;
  }
  cout << endl;
}
template<typename T, class Func>
void Solver<T, Func>::swap_rows(Matrix<T>& A, int row1, int row2) {
  // cout << "Swapped rows " << row1 << " " << row2 << endl;
  if (row1 == row2) return;
  swap(A[row1], A[row2]);
}
template<typename T, class Func>
void Solver<T, Func>::swap_columns(Matrix<T>& A, int col1, int col2, vector<int>& ans_order) {
  // cout << "Swapped columns " << col1 << " " << col2 << endl;
  swap(ans_order[col1], ans_order[col2]);
  if (col1 == col2) return;
  for (int row = 0; row < A.size(); ++row) {
    swap(A[row][col1], A[row][col2]);
  }
}
template<typename T, class Func>
void Solver<T, Func>::swap_columns(Matrix<T>& A, int col1, int col2) {
  // cout << "Swapped columns " << col1 << " " << col2 << endl;
  if (col1 == col2) return;
  for (int row = 0; row < A.size(); ++row) {
    swap(A[row][col1], A[row][col2]);
  }
}
template<typename T, class Func>
void Solver<T, Func>::normalize_row(Matrix<T>& A, int row) {
  T koef = A[row][row];
  for (int col = row; col < A[0].size(); ++col) {
    A[row][col] /= koef;
  }
}
template<typename T, class Func>
Matrix<T> Solver<T, Func>::SolveSystem(Matrix<T> A, Matrix<T> B, const SolverMethod& method_) {
  auto AA(A);
  vector<int> ans_order(A.size());
  for (int i = 0; i < A.size(); ++i) {
    ans_order[i] = i;
  }

  A = makeExtended(A, B);

  // Print(A, ans_order);

  for (int stage = 0; stage < A.size(); ++stage) {
    SolveStage(A, stage, method_, ans_order);
    // Print(A, ans_order);
  }

  // Print(A, ans_order);

  for (int i = A.size() - 1; i >= 0; --i) {
    for (int j = i + 1; j < A.size(); ++j) {
      // B[i] -= A[i][j] * B[j];
      sub_row(A, i, j, A[i][j]);
      A[i][j] = T(0);
    }
  }

  // Print(A, ans_order);


  auto ans = extractAnsMatrixFromExtended(A);
  Matrix<T> X(vector<vector<T >>(ans.size(), vector<T>(ans[0].size(), T(0))));
  for (int i = 0; i < ans.size(); ++i) {
    X[ans_order[i]] = ans[i];
  }
  // assert(mult(AA, X) == B);
  return ans;
}

template<typename T, class Func>
void Solver<T, Func>::SolveStage(Matrix<T>& A, int stage, const SolverMethod& method_, vector<int>& ans_order) {
  // cout << "Solve stage(for A) " << stage << endl;
  switch (method_) {
    case SolverMethod::DO_NOT_TOUCH : {
      for (int row = stage; row < A.size(); ++row) {
        if (A[row][stage] != T(0)) {
          swap_rows(A, stage, row);
          break;
        }
      }
      break;
    }
    case SolverMethod::BEST_IN_COLUMN : {
      swap_rows(A, stage, best_in_the_col(A, stage));
      break;
    }
    case SolverMethod::BEST_IN_ROW : {
      swap_columns(A, stage, best_in_the_row(A, stage), ans_order);
      break;
    }
    case SolverMethod::BEST_IN_MATRIX : {
      auto max_pos = best_in_the_sqr(A, stage, stage);
      swap_rows(A, stage, max_pos.first);
      swap_columns(A, stage, max_pos.second, ans_order);
      break;
    }
  }
  if (A[stage][stage] == T(0)) return;

  normalize_row(A, stage);
  for (int row = stage + 1; row < A.size(); ++row) {
    substract_str(A, row, stage);
  }

}
template<typename T, class Func>
Matrix<T> Solver<T, Func>::makeExtended(const Matrix<T>& A, const Matrix<T>& B) {
  Matrix<T> AA(A);
  for (int i = 0; i < AA.size(); ++i) {
    std::copy(begin(B[i]), end(B[i]), std::back_inserter(AA[i]));
  }
  return AA;
}
template<typename T, class Func>
void Solver<T, Func>::sub_row(Matrix<T>& A, int row1, int row2, T koef) {
  for (int j = 0; j < A[0].size(); ++j) {
    A[row1][j] -= A[row2][j] * koef;
  }
}
template<typename T, class Func>
Matrix<T> Solver<T, Func>::extractAnsMatrixFromExtended(const Matrix<T>& A) {
  Matrix<T> B(A.size(), vector<T>(A[0].size() - A.size()));
  for (int i = 0; i < A.size(); ++i) {
    for (int j = A.size(); j < A[0].size(); ++j) {
      B[i][j - A.size()] = A[i][j];
    }
  }
  return B;
}
template<typename T, class Func>
tuple<Matrix<T>, Matrix<T>, Matrix<T>, Matrix<T>> Solver<T, Func>::DLUP_Step(Matrix<T>& A, int stage) {

  Matrix<T> I(getIdentityMatrix<T>(A.size()));
  Matrix<T> P(I), D(I), L(I);
  auto best_pos = best_in_the_sqr(A, stage, stage);

  swap(D[stage], D[best_pos.first]);
  swap(P[stage], P[best_pos.second]);

  swap_rows(A, stage, best_pos.first);
  vector<int> ans_order(A.size());
  swap_columns(A, stage, best_pos.second, ans_order);

  L[stage][stage] = T(1) / A[stage][stage];
  for (int col = stage; col < A.size(); ++col) {
    A[stage][col] /= A[stage][stage];
  }

  for (int row = stage + 1; row < A.size(); ++row) {
    L[row][stage] = A[row][stage] / A[stage][stage];
    sub_row(A, row, stage, A[row][stage] / A[stage][stage]);
  }

  return tie(L, D, A, P);
}
template<typename T, class Func>
tuple<Matrix<T>, Matrix<T>, Matrix<T>, Matrix<T>> Solver<T, Func>::DLUP_Decomposition(Matrix<T> A,
                                                                                      const SolverMethod& method_) {
  int n = A.size();
  vector<Matrix<T>> Ps, Ds;
  Matrix<T> L(n, vector<T>(n, T(0))), U(move(A)), I(getIdentityMatrix<T>(n));

  // vector<pair<int, int>>

  bool is_first = true;
  for (int stage = 0; stage < n; ++stage) {
    pair<int, int> best_pos;
    switch (method_) {
      case SolverMethod::BEST_IN_MATRIX: {
        best_pos = best_in_the_sqr(U, stage, stage);
        break;
      }
      case SolverMethod::BEST_IN_ROW: {
        best_pos = {stage, best_in_the_row(U, stage)};
        break;
      }
      case SolverMethod::BEST_IN_COLUMN: {
        best_pos = {best_in_the_col(U, stage), stage};
        break;
      }
      case SolverMethod::DO_NOT_TOUCH: {
        best_pos = {stage, stage};
        break;
      }
      default:throw logic_error("TBD");
    }

    swap_rows(U, stage, best_pos.first);
    swap_columns(U, stage, best_pos.second);

//    if (stage != best_pos.first) {
//      cout << "change rows fact " << stage << " " << best_pos.first << "\n";
//    }
//    if (stage != best_pos.second) {
//      cout << "change column fact " << stage << " " << best_pos.second << "\n";
//    }

    Matrix<T> P(I), D(I);
    swap(P[stage], P[best_pos.second]);
    swap(D[stage], D[best_pos.first]);
    Ps.push_back(P);
    Ds.push_back(D);

    // PrintMatrix(U, "U after changes");

    if (is_first) {
      is_first = false;
    } else {
      L = mult(D, L);
    }

    for (int row = stage; row < n; ++row) {
      L[row][stage] = U[row][stage];
    }

    assert(U[stage][stage] != T(0));
    T k = U[stage][stage];
    for (int col = stage; col < n; ++col) {
      U[stage][col] /= k;
    }

    for (int row = stage + 1; row < n; ++row) {
      sub_row(U, row, stage, U[row][stage]);
    }

//    cout << "After stage " << stage << endl;
//    PrintMatrix(U, "U" + to_string(stage));
//    PrintMatrix(L, "L" + to_string(stage));
//    PrintMatrix(P, "P" + to_string(stage));
//    PrintMatrix(D, "D" + to_string(stage));
  }

  Matrix<T> P = I;
  for (int i = 0; i < Ps.size(); ++i) {
    P = mult(P, Ps[i]);
  }

  Matrix<T> D = I;
  for (int i = Ds.size() - 1; i >= 0; --i) {
    D = mult(D, Ds[i]);
  }

  // D = GetReversed(D);
  // P = GetReversed(P);

  return tie(D, L, U, P);
}

template<typename T, class Func>
Matrix<T> Solver<T, Func>::GetReversed(Matrix<T> A) {
  auto ans = SolveSystem(A, getIdentityMatrix<T>(A.size()), SolverMethod::DO_NOT_TOUCH);
  // assert(mult(A, ans) == getIdentityMatrix<T>(A.size()));
  return ans;
}

template<typename T, class Func>
tuple<Matrix<T>, vector<T>> Solver<T, Func>::LDL_Decomposition(Matrix<T> A) {
  vector<T> D(A.size(), T(1));
  for (int now = 0; now < A.size(); ++now) {
    if (A[now][now] == 0) {
      throw std::invalid_argument("bad matrix");
    }
    T k;
    if (A[now][now] < T(0)) {
      k = sqrt(-A[now][now]);
      D[now] = T(-1);
    } else {
      k = sqrt(A[now][now]);
    }
    for (int col = now; col < A.size(); ++col) {
      A[now][col] /= k;
    }
    for (int row = now + 1; row < A.size(); ++row) {
      sub_row(A, row, now, A[row][now] / A[now][now]);
    }

  }
  A = Transpose(A);
  return tie(A, D);
}

template<typename T, class Func>
vector<T> Solver<T, Func>::SolveThreeDiagonalSystem(Matrix<T> A, vector<T> b) {
  for (int i = 1; i < A.size(); ++i) {
    if (A[i - 1][1] == 0) throw std::logic_error("divide by zero");
    T k = A[i][0] / A[i - 1][1];

    A[i][0] -= k * A[i - 1][1];
    A[i][1] -= k * A[i - 1][2];
    b[i] -= k * b[i - 1];
  }

  // PrintThreeDiagonal(A);
  // PrintMatrix(b, "b");

  b.back() /= A.back()[1];
  A.back()[1] = 1;
  for (int i = A.size() - 2; i >= 0; --i) {
    b[i] -= A[i][2] * b[i + 1];
    b[i] /= A[i][1];
    A[i][1] = 1;
  }
  PrintThreeDiagonal(A);
  return b;
}
template<typename T, class Func>
void Solver<T, Func>::PrintThreeDiagonal(const Matrix<T>& A) {
  int n = A.size();
  cout << A[0][1] << " " << A[0][2] << " ";
  for (int i = 0; i < n - 2; ++i) {
    cout << T(0) << " ";
  }
  cout << endl;
  int k = 0;
  for (int i = 1; i < n - 1; ++i) {
    for (int j = 0; j < k; ++j) {
      cout << T(0) << " ";
    }
    cout << A[i][0] << " ";
    cout << A[i][1] << " ";
    cout << A[i][2] << " ";
    for (int j = 0; j < n - 3 - k; ++j) {
      cout << T(0) << " ";
    }
    cout << endl;
    ++k;
  }
  for (int j = 0; j < n - 2; ++j) {
    cout << T(0) << " ";
  }
  cout << A.back()[0] << " " << A.back()[1] << endl;
}

template<typename T, class Func>
vector<T> Solver<T, Func>::SolveUp(const Matrix<T>& A, const vector<T>& b) {
  assert(A.size() == A[0].size() && b.size() == A.size());
  int n = A.size();
  vector<T> ans(b);
  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j) {
      ans[i] -= ans[j] * A[i][j];
    }
    ans[i] /= A[i][i];
  }
  return ans;
}

template<typename T, class Func>
vector<T> Solver<T, Func>::SolveDown(const Matrix<T>& A, const vector<T>& b) {
  assert(A.size() == A[0].size() && b.size() == A.size());
  int n = A.size();
  vector<T> ans(b);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      ans[i] -= ans[j] * A[i][j];
    }
    ans[i] /= A[i][i];
  }
  return ans;
}

template<typename T, class Func>
vector<T> Solver<T, Func>::SolveLinearSystemUsingLDLt(Matrix<T> A, vector<T> B) {

  Matrix<T> L;
  vector<T> D;
  tie(L, D) = LDL_Decomposition(A);
  // L * D * L^T == A

  // PrintMatrix(D, "D");
  // PrintMatrix(L, "L");
  // PrintMatrix(Transpose(L), "L^T");
  // PrintMatrix(mult(L, Transpose(L)), "L /cdot L^T");
  // Matrix<T> DD(A.size(), vector<T>(A.size(), 0));
  // for (int i = 0; i < A.size(); ++i) {
  //   DD[i][i] = D[i];
  // }
  // PrintMatrix(mult(mult(L, DD), Transpose(L)), "LDLt");

  // L D Lt X = B
  // L y = b
  // (D Lt) x = y

  int n = A.size();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      B[i] -= L[i][j] * D[j] * B[j];
    }
    B[i] /= L[i][i] * D[i];
  }

  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j) {
      B[i] -= L[j][i] * B[j];
    }
    B[i] /= L[i][i];
  }
  return B;
}

template<typename T, class Func>
Matrix<T> Solver<T, Func>::Transpose(Matrix<T> A) {
  assert(A.size() == A[0].size());
  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      swap(A[i][j], A[j][i]);
    }
  }
  return A;
}
template<typename T, class Func>
vector<T> Solver<T, Func>::SolveLinearSystemUsingDLUP(Matrix<T> A, vector<T> b) {
  int n = A.size();

  Matrix<T> D, L, U, P;
  tie(D, L, U, P) = DLUP_Decomposition(A, SolverMethod::BEST_IN_COLUMN);
  PrintMatrix(D, "D");
  PrintMatrix(L, "L");
  PrintMatrix(U, "U");
  PrintMatrix(P, "P");
  if (!compareMatrixOfDouble(A, mult(mult(GetReversed(D), L), mult(U, GetReversed(P))))) {
    cerr << "Smth goes wrong" << endl;
    assert(false);
  }
  vector<T> b_original(b);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (D[i][j]) {
        b[i] = b_original[j];
      }
    }
  }
  // L Y = B
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      b[i] -= L[i][j] * b[j];
    }
    b[i] /= L[i][i];
  }

  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j) {
      b[i] -= U[i][j] * b[j];
    }
  }

  return b;
}
