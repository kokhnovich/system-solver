
#include "MatrixChanger.h"
template<typename T>
bool checkUpperMatrix(const Matrix<T>& a) {
  assert(a.size() == a[0].size());
  for (int i = 0; i < a.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      assert(a[i][j] == T(0));
    }
  }
}
template<typename T>
bool checkLowerMatrix(const Matrix<T>& a) {
  assert(a.size() == a[0].size());
  for (int i = 0; i < a.size(); ++i) {
    for (int j = i + 1; j < a.size(); ++j) {
      assert(a[i][j] == T(0));
    }
  }
}

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

}

void SolverDoubleTests() {

}

void SolverFractionTests() {

}

void TestLU() {

}

void TestMatrixChanger() {
  Matrix<double> m = {{1, 2, 3},
                      {4, 4, 5},
                      {5, 6, 2}};
  MatrixChanger matrix_changer(m);
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
}
void RunTests() {
  FractionTests();
}