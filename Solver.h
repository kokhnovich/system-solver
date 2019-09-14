//
// Created by user on 13.09.19.
//

#ifndef SYSTEM_SOLVER__SOLVER_H_
#define SYSTEM_SOLVER__SOLVER_H_

enum class SolverMethod {
  DO_NOT_TOUCH,
  BEST_IN_ROW,
  BEST_IN_COLUMN,
  BEST_IN_MATRIX
};

template<typename T>
vector<T>& operator-=(vector<T>& one, const vector<T>& two) {
  if (one.size() != two.size()) throw logic_error("a -= b bad sizes");
  for (size_t i = 0; i < one.size(); ++i) {
    one[i] -= two[i];
  }
  return one;
}

template<typename T>
vector<T> operator*(T koef, const vector<T>& one) {
  vector<T> ans(one);
  for (auto& i : ans) {
    i = i * koef;
  }
  return ans;
}

/** matrix should be NxN
 *  template class Func can be, for example, std:greater<T> or std:less<T>
 *  Func must be linear!
 *  T for elements of A
 *
 *  @TODO
 *  //T2 for b: T2=T , if b is a vector
 *  //          T2=vector<T> if b is a matrix // you should also implement operator-= for two vector<T>
 *
**/
template<typename T, class Func=std::greater<T>>
class Solver {
 private:

 public:
  explicit Solver() = default;

  tuple<vector<vector<T>>, vector<vector<T>>, vector<vector<T>>> LUP_Decomposition(const vector<vector<T>>& A,
                                                                                   const SolverMethod& method_);
  // vector<T> SolveSystemUsingLU();
  vector<vector<T>> SolveSystem(vector<vector<T>> A, vector<vector<T>> B, const SolverMethod& method_);

  void Print(const vector<vector<T>>& A, const vector<int>& ans_order) const;

 private:

  void SolveStage(vector<vector<T>>& A, int stage, const SolverMethod& method_, vector<int>& ans_order);

  vector<vector<T>> makeExtended(const vector<vector<T>>& A, const vector<vector<T>>& B);
  vector<vector<T>> extractAnsMatrixFromExtended(const vector<vector<T>>& A);

  int best_in_the_row(const vector<vector<T>>& A, int row);
  int best_in_the_col(const vector<vector<T>>& A, int col);
  pair<int, int> best_in_the_sqr(const vector<vector<T>>& A, int start_i, int start_j);

  void swap_rows(vector<vector<T>>& A, int row1, int row2);
  void swap_columns(vector<vector<T>>& A, int col1, int col2, vector<int>& ans_order);

  void normalize_row(vector<vector<T>>& A, int row);
  void substract_str(vector<vector<T>>& A, int row, int stage);
  void sub_row(vector<vector<T>>& A, int row1, int row2, T koef); // row1 -= row2 * koef
};

template<typename T, class Func>
tuple<vector<vector<T>>, vector<vector<T>>, vector<vector<T>>> Solver<T, Func>::LUP_Decomposition(const vector<vector<T>>& A,
                                                                                                  const SolverMethod& method_) {
  vector<int> ans_order(A.size());
  vector<vector<T>> L(A.size(), vector<T>(A.size(), 0)), U(A), P(A.size(), vector<T>(A.size(), 0));
  for (int i = 0; i < A.size(); ++i) {
    ans_order[i] = i;
    P[i][i] = 1;
  }
  Print(U, ans_order);
  for (int stage = 0; stage < U.size(); ++stage) {
    L[stage][stage] = 1;
    switch (method_) {
      case SolverMethod::BEST_IN_ROW : {
        swap_columns(U, stage, best_in_the_row(U, stage), ans_order);
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
      // Print(A, ans_order);
    }

    Print(U, ans_order);
  }

//  for (int i = 0; i < A.size(); ++i) {
//
//  }
//  vector<vector<T>> ans(size_, vector<T>(size_, T(0)));
//  for (int i = 0; i < size_; ++i) {
//    P[i][ans_order[i]] = 1;
//      for (int j = 0; j < size_; ++j) {
//        ans[j][ans_order[i]] = a[j][i];
//      }
//  }
  return make_tuple(L, U, P);
}

template<typename T, class Func>
void Solver<T, Func>::substract_str(vector<vector<T>>& A, int row, int stage) {
  T koef = A[row][stage];
  for (int col = stage; col < A[0].size(); ++col) {
    A[row][col] -= A[stage][col] * koef;
  }
}

/**
template<typename T, class Func>
vector<T> Solver<T, Func>::SolveSystemUsingLU() {
  if (!was_lu) throw logic_error("You should do LU before that");
//
//   * A = LU
//   * Ax=b
//   * LUx=b
//   * [y=Ux]
//   * 1) Ly=b, L is lower-diagonal
//   * 2) Ux=y, L is upper-diagonal
//   * return x
//
  vector<T> y(size_, T(0));

  // a = mult(a, P);
  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < size_; ++j) {
      a[j][ans_order[i]] = a[j][i];
    }
  }

  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < i; ++j) {
      b[i] -= y[j] * L[i][j];
    }
    y[i] = b[i];
  }

  vector<T> x(size_, T(0));
  for (int i = size_ - 1; i >= 0; --i) {
    for (int j = i + 1; j < size_; ++j) {
      y[i] -= a[i][j] * x[j];
    }
    x[i] = y[i] / a[i][i];
  }

  return x;
}
**/

template<typename T, class Func>
int Solver<T, Func>::best_in_the_row(const vector<vector<T>>& A, int row) {
  int best_column = row;
  for (int column = row + 1; column < A[0].size(); ++column) {
    if (Func()(A[row][column], A[row][best_column])) {
      best_column = column;
    }
  }
  return best_column;
}
template<typename T, class Func>
int Solver<T, Func>::best_in_the_col(const vector<vector<T>>& A, int col) {
  int best_row = col;
  for (int i = col + 1; i < A[0].size(); ++i) {
    if (Func()(A[i][col], A[best_row][col])) {
      best_row = i;
    }
  }
  return best_row;
}
template<typename T, class Func>
pair<int, int> Solver<T, Func>::best_in_the_sqr(const vector<vector<T>>& A, int start_i, int start_j) {
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
void Solver<T, Func>::Print(const vector<vector<T>>& A, const vector<int>& ans_order) const {
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
void Solver<T, Func>::swap_rows(vector<vector<T>>& A, int row1, int row2) {
  cout << "Swapped rows " << row1 << " " << row2 << endl;
  if (row1 == row2) return;
  swap(A[row1], A[row2]);
}
template<typename T, class Func>
void Solver<T, Func>::swap_columns(vector<vector<T>>& A, int col1, int col2, vector<int>& ans_order) {
  cout << "Swapped columns " << col1 << " " << col2 << endl;
  swap(ans_order[col1], ans_order[col2]);
  if (col1 == col2) return;
  for (int row = 0; row < A.size(); ++row) {
    swap(A[row][col1], A[row][col2]);
  }
}
template<typename T, class Func>
void Solver<T, Func>::normalize_row(vector<vector<T>>& A, int row) {
  T koef = A[row][row];
  for (int col = row; col < A[0].size(); ++col) {
    A[row][col] /= koef;
  }
}
template<typename T, class Func>
vector<vector<T>> Solver<T, Func>::SolveSystem(vector<vector<T>> A, vector<vector<T>> B, const SolverMethod& method_) {

  vector<int> ans_order(A.size());
  for (int i = 0; i < A.size(); ++i) {
    ans_order[i] = i;
  }

  A = makeExtended(A, B);

  Print(A, ans_order);

  for (int stage = 0; stage < A.size(); ++stage) {
    SolveStage(A, stage, method_, ans_order);
    Print(A, ans_order);
  }

  Print(A, ans_order);

  for (int i = A.size() - 1; i >= 0; --i) {
    for (int j = i + 1; j < A.size(); ++j) {
      // B[i] -= A[i][j] * B[j];
      sub_row(A, i, j, A[i][j]);
      A[i][j] = T(0);
    }
  }

  Print(A, ans_order);

//  vector<vector<T>> ans(A.size());
//  for (size_t i = 0; i < A.size(); ++i) {
//    ans[ans_order[i]] = B[i];
//  }
  return extractAnsMatrixFromExtended(A);
}

template<typename T, class Func>
void Solver<T, Func>::SolveStage(vector<vector<T>>& A, int stage, const SolverMethod& method_, vector<int>& ans_order) {
  cout << "Solve stage(for A) " << stage << endl;
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
vector<vector<T>> Solver<T, Func>::makeExtended(const vector<vector<T>>& A, const vector<vector<T>>& B) {
  vector<vector<T>> AA(A);
  for (int i = 0; i < AA.size(); ++i) {
    std::copy(begin(B[i]), end(B[i]), std::back_inserter(AA[i]));
  }
  return AA;
}
template<typename T, class Func>
void Solver<T, Func>::sub_row(vector<vector<T>>& A, int row1, int row2, T koef) {
  for (int j = 0; j < A[0].size(); ++j) {
    A[row1][j] -= A[row2][j] * koef;
  }
}
template<typename T, class Func>
vector<vector<T>> Solver<T, Func>::extractAnsMatrixFromExtended(const vector<vector<T>>& A) {
  vector<vector<T>> B(A[0].size() - A.size(), vector<T>(A[0].size() - A.size()));
  for (int i = 0; i < A.size(); ++i) {
    for (int j = A.size(); j < A[0].size(); ++j) {
      B[i][j - A.size()] = A[i][j];
    }
  }
  return B;
}

#endif //SYSTEM_SOLVER__SOLVER_H_
