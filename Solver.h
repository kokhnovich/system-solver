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
  bool was_lu;
  vector<vector<T>> a, L, P;
  vector<T> b;
  vector<int> ans_perm;
  int size_;
  SolverMethod method_;

 public:
  explicit Solver(vector<vector<T>> aa, vector<T> bb, const SolverMethod& method = SolverMethod::DO_NOT_TOUCH)
      : a(move(aa)), b(move(bb)), was_lu(false), ans_perm(vector<int>(a.size(), 0)),
        method_(method), L(vector<vector<T>>(a.size(), vector<T>(a.size(), T(0)))),
        P(vector<vector<T>>(a.size(), vector<T>(a.size(), T(0)))) {
    if (a.empty() || a.size() != a[0].size() || a.size() != b.size()) throw logic_error("bad matrix");
    size_ = a.size();
    for (int i = 0; i < size_; ++i) {
      ans_perm[i] = i;
    }
  }

  tuple<vector<vector<T>>, vector<vector<T>>, vector<vector<T>>> LUP_Decomposition();
  vector<T> SolveSystemUsingLU();
  vector<T> GetSolution();
  vector<T> SolveSystem(vector<vector<T>> A, vector<T> B);

  void Print() const;

 private:

  void SolveStage(int stage, bool lu = false);

  int best_in_the_row(int row);
  int best_in_the_col(int col);
  pair<int, int> best_in_the_sqr(int start_i, int start_j);

  void swap_rows(int row1, int row2);
  void swap_columns(int col1, int col2);

  void normalize_row(int row, bool work_with_ans = true);
  void substract_str(int row, int stage, bool work_with_ans = true);
};

template<typename T, class Func>
vector<T> Solver<T, Func>::GetSolution() {
  /// time complexity is O(n^3)
  for (int i = 0; i < size_; ++i) {
    SolveStage(i);
    cout << "stage " << i << endl;
    Print();
  }

  //Print();
  //cout << endl;

  /// time complexity is O(n^2)
  for (int i = size_ - 1; i >= 0; --i) {
    T cost = T(0);
    for (int j = i + 1; j < size_; ++j) {
      cost += a[i][j] * b[j];
      a[i][j] = T(0);
    }
    b[i] -= cost;
  }
  vector<T> ans(size_);
  for (int i = 0; i < size_; ++i) {
    ans[ans_perm[i]] = b[i];
  }
  return ans;
}

template<typename T, class Func>
tuple<vector<vector<T>>, vector<vector<T>>, vector<vector<T>>> Solver<T, Func>::LUP_Decomposition() {
  for (int i = 0; i < size_; ++i) {
    SolveStage(i, true);
    // cout << "stage " << i << endl;
    // Print();
  }
  Print();
  vector<vector<T>> ans(size_, vector<T>(size_, T(0)));
  for (int i = 0; i < size_; ++i) {
    P[i][ans_perm[i]] = 1;
//      for (int j = 0; j < size_; ++j) {
//        ans[j][ans_perm[i]] = a[j][i];
//      }
  }

  was_lu = true;
  return make_tuple(L, a, P);
}

template<typename T, class Func>
void Solver<T, Func>::substract_str(int row, int stage, bool work_with_ans) {
  T koef;
  if (work_with_ans) {
    koef = a[row][stage];
  } else {
    koef = a[row][stage] / a[stage][stage];
    int diff = row - stage;
    L[row][stage] = koef;
  }
  for (int col = stage; col < size_; ++col) {
    a[row][col] -= a[stage][col] * koef;
  }
  if (work_with_ans) {
    b[row] -= b[stage] * koef;
  }
}
template<typename T, class Func>
void Solver<T, Func>::SolveStage(int stage, bool lu) {
  cout << "Solve stage " << stage << endl;
  switch (method_) {
    case SolverMethod::DO_NOT_TOUCH : {
      if (lu) break;
      for (int row = stage; row < size_; ++row) {
        if (a[row][stage] != T(0)) {
          swap_rows(stage, row);
          break;
        }
      }
      break;
    }
    case SolverMethod::BEST_IN_COLUMN : {
      swap_rows(stage, best_in_the_col(stage));
      break;
    }
    case SolverMethod::BEST_IN_ROW : {
      swap_columns(stage, best_in_the_row(stage));
      break;
    }
    case SolverMethod::BEST_IN_MATRIX : {
      auto max_pos = best_in_the_sqr(stage, stage);
      swap_rows(stage, max_pos.first);
      swap_columns(stage, max_pos.second);
      break;
    }
  }
  if (a[stage][stage] == T(0)) return;

  if (!lu) {
    normalize_row(stage);
    for (int row = stage + 1; row < size_; ++row) {
      substract_str(row, stage);
    }
  } else {
    L[stage][stage] = 1;
    for (int row = stage + 1; row < size_; ++row) {
      substract_str(row, stage, false);
    }
  }
}

template<typename T, class Func>
vector<T> Solver<T, Func>::SolveSystemUsingLU() {
  if (!was_lu) throw logic_error("You should do LU before that");
  /*
   * A = LU
   * Ax=b
   * LUx=b
   * [y=Ux]
   * 1) Ly=b, L is lower-diagonal
   * 2) Ux=y, L is upper-diagonal
   * return x
   */

  vector<T> y(size_, T(0));

  // a = mult(a, P);
  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < size_; ++j) {
      a[j][ans_perm[i]] = a[j][i];
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
template<typename T, class Func>
int Solver<T, Func>::best_in_the_row(int row) {
  int best_column = row;
  for (int column = row + 1; column < size_; ++column) {
    if (Func()(a[row][column], a[row][best_column])) {
      best_column = column;
    }
  }
  return best_column;
}
template<typename T, class Func>
int Solver<T, Func>::best_in_the_col(int col) {
  int best_row = col;
  for (int i = col + 1; i < size_; ++i) {
    if (Func()(a[i][col], a[best_row][col])) {
      best_row = i;
    }
  }
  return best_row;
}
template<typename T, class Func>
pair<int, int> Solver<T, Func>::best_in_the_sqr(int start_i, int start_j) {
  int best_i = start_i, best_j = start_j;
  for (int i = start_i; i < size_; ++i) {
    for (int j = start_j; j < size_; ++j) {
      if (Func()(a[i][j], a[best_i][best_j])) {
        best_i = i, best_j = j;
      }
    }
  }
  return make_pair(best_i, best_j);
}
template<typename T, class Func>
void Solver<T, Func>::Print() const {
  for (int i = 0; i < size_; ++i) {
    for (auto& j : a[i]) {
      cout << j << " ";
    }
    cout << b[i] << endl;
  }
  cout << "Порядок неизвестных: ";
  for (const auto& i : ans_perm) {
    cout << i + 1 << " ";
  }
  cout << endl;
}
template<typename T, class Func>
void Solver<T, Func>::swap_rows(int row1, int row2) {
  cout << "Swapped rows " << row1 << " " << row2 << endl;
  if (row1 == row2) return;
  swap(a[row1], a[row2]);
  swap(b[row1], b[row2]);
}
template<typename T, class Func>
void Solver<T, Func>::swap_columns(int col1, int col2) {
  cout << "Swapped columns " << col1 << " " << col2 << endl;
  if (col1 == col2) return;
  swap(ans_perm[col1], ans_perm[col2]);
  for (int row = 0; row < size_; ++row) {
    swap(a[row][col1], a[row][col2]);
  }
}
template<typename T, class Func>
void Solver<T, Func>::normalize_row(int row, bool work_with_ans) {
  T koef = a[row][row];
  for (int col = row; col < size_; ++col) {
    a[row][col] /= koef;
  }
  if (work_with_ans) {
    b[row] /= koef;
  }
}
template<typename T, class Func>
vector<T> Solver<T, Func>::SolveSystem(vector<vector<T>> A, vector<T> B) {


}


#endif //SYSTEM_SOLVER__SOLVER_H_
