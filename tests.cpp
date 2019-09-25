
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

void RunTests() {
  FractionTests();
}