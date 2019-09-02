#include <bits/stdc++.h>

#include <utility>

using namespace std;

class Fraction {
 private:
  int numerator, denominator;

 public:
  // @TODO negatives
  Fraction(int n, int d) : numerator(n), denominator(d) {
    if (denominator == 0) throw logic_error("divide by zero");
    Normalize();
  }

  void Normalize() {
    if (numerator * denominator < 0) {
      numerator = -1 * abs(numerator);
      denominator = abs(denominator);
    } else {
      numerator = abs(numerator);
      denominator = abs(denominator);
    }
    int gcd = __gcd(abs(numerator), denominator);
    numerator /= gcd;
    denominator /= gcd;
  }

  Fraction operator+(const Fraction& otherFraction) const {
    int n = numerator * otherFraction.denominator + otherFraction.numerator * denominator;
    int d = denominator * otherFraction.denominator;
    return {n, d};
  }
  Fraction operator-(const Fraction& otherFraction) const {
    int n = numerator * otherFraction.denominator - otherFraction.numerator * denominator;
    int d = denominator * otherFraction.denominator;
    return {n, d};
  }
  Fraction operator*(const Fraction& otherFraction) const {
    int n = numerator * otherFraction.numerator;
    int d = denominator * otherFraction.denominator;
    return {n, d};
  }
  Fraction operator/(Fraction& otherFraction) const {
    int n = numerator * otherFraction.denominator;
    int d = denominator * otherFraction.numerator;
    return {n, d};
  }

  bool operator==(const Fraction& otherFraction) const {
    return numerator == otherFraction.numerator && denominator == otherFraction.denominator;
  }
  bool operator!=(const Fraction& otherFraction) const {
    return !(*this == otherFraction);
  }

  void show() {
    cout << numerator << "/" << denominator << endl;
  }
};

template<typename T>
class Solver {
 public:
  /// matrix should be NxN
  // @TODO add ability to change function, default=std::max
  // @TODO add METHOD param
  // @TODO zero assertion
  explicit Solver(vector<vector<int>> aa, vector<int> bb)
      : a(move(aa)), b(move(bb)), ans_perm(vector<int>(a.size(), 0)) {
    if (a.empty() || a.size() != a[0].size()) throw logic_error("bad matrix error");
    size_ = a.size();
    for (int i = 0; i < size_; ++i) {
      ans_perm[i] = i;
    }
  }

  void Solve() {
    // up-down
    for (int i = 0; i < size_; ++i) {
      Solve(i, i);
    }
    Print();
    cout << endl;

    // down-up
    for (int i = size_ - 1; i >= 0; --i) {
      int cost = 0;
      for (int j = i + 1; j < size_; ++j) {
        cost += a[i][j] * b[j];
        a[i][j] = 0;
      }
      b[i] -= cost;
    }
  }

  void Solve(int start_i, int start_j) {
    //auto max_pos = max_in_the_sqr(start_i, start_j);
    auto max_pos = make_pair(start_i, start_j);
    //swap_rows(start_i, max_pos.first);
    //swap_columns(start_j, max_pos.second);
    normalize_row(start_i, start_i);
    for (int row = start_i + 1; row < size_; ++row) {
      substract_str(row, start_i, start_i);
    }

  }

  /// returns number of column, where max is located
  const int max_in_the_row(int row) {
    int best_column = 0;
    for (int column = 1; column < size_; ++column) {
      if (a[row][column] > a[row][best_column]) {
        best_column = column;
      }
    }
    return best_column;
  }

  const pair<int, int> max_in_the_sqr(int start_i, int start_j) {
    int best_i = start_i, best_j = start_j;
    for (int i = start_i; i < size_; ++i) {
      for (int j = start_j; j < size_; ++j) {
        if (a[i][j] > a[best_i][best_j]) {
          best_i = i, best_j = j;
        }
      }
    }
    return make_pair(best_i, best_j);
  }

  void swap_rows(int row1, int row2) {
    swap(a[row1], a[row2]);
    swap(b[row1], b[row2]);
  }

  void swap_columns(int col1, int col2) {
    swap(ans_perm[col1], ans_perm[col2]);
    for (int row = 0; row < size_; ++row) {
      swap(a[row][col1], a[row][col2]);
    }
  }

  /// think that from-row is normalized
  void substract_str(int row, int from, int first) {
    int koef = a[row][first];
    for (int col = first; col < size_; ++col) {
      a[row][col] -= a[from][col] * koef;
    }
    b[row] -= b[from] * koef;
  }

  void normalize_row(int row, int first) {
    int koef = a[row][first];
    a[row][first] = 1;
    for (int col = 1; col < size_; ++col) {
      a[row][col] /= koef;
    }
  }

  void Print() {
    for (int i = 0; i < size_; ++i) {
      for (auto& j : a[i]) {
        cout << j << " ";
      }
      cout << b[i] << endl;
    }
    cout << endl;
  }

 private:
  vector<vector<int>> a;
  vector<int> b;
  vector<int> ans_perm;
  int size_;
};

void FractionTests() {
  Fraction a(1, 2), b(1, 3);
  assert(a + b == Fraction(5, 6));
  assert(b + a == Fraction(5, 6));
  assert(a - b == Fraction(1, 6));
  assert(b - a == Fraction(-1, 6));
  assert(b - a == Fraction(1, -6));
  assert(b - a != Fraction(1, 6));

  assert(a * b == Fraction(1, 6));
  assert(a / b == Fraction(3, 2));
  assert(a / b == Fraction(6, 4));
  cout << "Fraction passed tests" << endl;
}

// @TODO unit-tests

int main() {
  FractionTests();
  int n = 5;
  vector<vector<int>> a = {{1, 2, 3}, {4, 9, 15}, {3, 5, 7}};
  vector<int> b = {4, 19, 10};
  auto smth = new Solver<int>(a, b);
  smth->Print();
  smth->Solve();
  smth->Print();
  return 0;
}