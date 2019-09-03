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

/** matrix should be NxN
 *  template class Func can be, for example, std:greater<T> or std:less<T>
 *
**/
template<typename T, class Func=std::greater<T>>
class Solver {
 public:
  // @TODO add ability to change function, default=std::max
  // @TODO add METHOD param
  // @TODO zero assertion
  explicit Solver(vector<vector<T>> aa, vector<T> bb, const SolverMethod& method = SolverMethod::DO_NOT_TOUCH)
      : a(move(aa)), b(move(bb)), ans_perm(vector<int>(a.size(), 0)), method_(method) {
    if (a.empty() || a.size() != a[0].size() || a.size() != b.size()) throw logic_error("bad matrix");
    size_ = a.size();
    for (int i = 0; i < size_; ++i) {
      ans_perm[i] = i;
    }
  }

  void Solve() {
    /// time complexity is O(n^3)
    for (int i = 0; i < size_; ++i) {
      Solve(i);
      //cout << "stage " << i << endl;
      //Print();
    }

    //Print();
    //cout << endl;

    /// time complexity is O(n^2)
    for (int i = size_ - 1; i >= 0; --i) {
      T cost = 0;
      for (int j = i + 1; j < size_; ++j) {
        cost += a[i][j] * b[j];
        a[i][j] = 0;
      }
      b[i] -= cost;
    }
  }

  void Solve(int stage) {
    switch (method_) {
      case SolverMethod::DO_NOT_TOUCH : {
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

    normalize_row(stage);
    for (int row = stage + 1; row < size_; ++row) {
      substract_str(row, stage);
    }

  }

  /// returns number of column, where max is located
  int best_in_the_row(int row) {
    int best_column = row;
    for (int column = row + 1; column < size_; ++column) {
      if (Func()(a[row][column], a[row][best_column])) {
        best_column = column;
      }
    }
    return best_column;
  }

  int best_in_the_col(int col) {
    int best_row = col;
    for (int i = col + 1; i < size_; ++i) {
      if (Func()(a[i][col], a[best_row][col])) {
        best_row = i;
      }
    }
    return best_row;
  }

  pair<int, int> best_in_the_sqr(int start_i, int start_j) {
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

  void swap_rows(int row1, int row2) {
    if (row1 == row2) return;
    swap(a[row1], a[row2]);
    swap(b[row1], b[row2]);
  }

  void swap_columns(int col1, int col2) {
    if (col1 == col2) return;
    swap(ans_perm[col1], ans_perm[col2]);
    for (int row = 0; row < size_; ++row) {
      swap(a[row][col1], a[row][col2]);
    }
  }

  void substract_str(int row, int stage) {
    T koef = a[row][stage];
    for (int col = stage; col < size_; ++col) {
      a[row][col] -= a[stage][col] * koef;
    }
    b[row] -= b[stage] * koef;
  }

  void normalize_row(int row) {
    T koef = a[row][row];
    for (int col = row; col < size_; ++col) {
      a[row][col] /= koef;
    }
    b[row] /= koef;
  }

  void Print() {
    for (int i = 0; i < size_; ++i) {
      for (auto& j : a[i]) {
        cout << j << " ";
      }
      cout << b[i] << endl;
    }
    cout << "ans_perm: ";
    for (const auto& i : ans_perm) {
      cout << i << " ";
    }
    cout << endl;
  }

 private:
  vector<vector<T>> a;
  vector<T> b;
  vector<int> ans_perm;
  int size_;
  SolverMethod method_;
};

void FractionTests() {
  Fraction a(1, 2), b(1, 3);
  assert(a + b == Fraction(5, 6));
  assert(b - a == Fraction(1, -6));
  assert(a * b == Fraction(1, 6));
  assert(a / b == Fraction(6, 4));
  assert(Fraction(3, 26) < Fraction(51, 52));
  cout << "Fraction passed tests" << endl;
}

void SolverIntTests() {
  vector<vector<int>> a = {{1, 2, 3}, {4, 9, 15}, {3, 5, 7}};
  vector<int> b = {4, 19, 10};
  auto smth = new Solver<int>(a, b);
  smth->Print();
  smth->Solve();
  smth->Print();
  cout << "Solver passed tests with int with no method" << endl;
  /*
  smth = new Solver<int>(a, b, SolverMethod::BEST_IN_COLUMN);
  smth->Print();
  smth->Solve();
  smth->Print();
  cout << "Solver passed tests with int with best_in_column" << endl;
   */
}

void SolverDoubleTests() {
  vector<vector<double>> a = {{1, 2, 3}, {4, 9, 15}, {3, 5, 7}};
  vector<double> b = {4, 19, 10};
  auto smth = new Solver<double>(a, b);
  smth->Print();
  smth->Solve();
  smth->Print();
  cout << "Solver passed tests with double" << endl;
}

void SolverFractionTests() {
  vector<vector<int>> a = {{1, 2, 3}, {4, 9, 15}, {3, 5, 7}};
  vector<int> b = {4, 19, 10};
  vector<vector<Fraction>> test_a(a.size(), vector<Fraction>(a.size()));
  vector<Fraction> test_b(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    for (size_t j = 0; j < a.size(); ++j) {
      test_a[i][j] = Fraction(a[i][j], 1);
    }
    test_b[i] = Fraction(b[i], 1);
  }

  auto solver = new Solver<Fraction>(test_a, test_b);
  solver->Print();
  solver->Solve();
  solver->Print();
  cout << "Solver passed tests with Fraction with no method" << endl;

  solver = new Solver<Fraction>(test_a, test_b, SolverMethod::BEST_IN_COLUMN);
  solver->Print();
  solver->Solve();
  solver->Print();
  cout << "Solver passed tests with Fraction with best_in_col" << endl;

  solver = new Solver<Fraction>(test_a, test_b, SolverMethod::BEST_IN_ROW);
  solver->Print();
  solver->Solve();
  solver->Print();
  cout << "Solver passed tests with Fraction with best_in_row" << endl;

  solver = new Solver<Fraction>(test_a, test_b, SolverMethod::BEST_IN_MATRIX);
  solver->Print();
  solver->Solve();
  solver->Print();
  cout << "Solver passed tests with Fraction with best_in_matrix" << endl;
}

int main() {
  FractionTests();
  //SolverIntTests();
  //SolverDoubleTests();
  SolverFractionTests();
  return 0;
}