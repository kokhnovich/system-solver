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

enum class DebugMethod {
  NO_DEBUG,
  MAIN_STEPS,
  ALL_STEPS
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

template<typename T>
struct ThreeDiagonal {
  T a, b, c;
};

template<typename T, class Func=greater_using_abs<T>>
class Solver {
 public:
  explicit Solver() = default;
  tuple<Matrix<T>, Matrix<T>, Matrix<T>> LUP_Decomposition(const Matrix<T>& A,
                                                           const SolverMethod& method_ = SolverMethod::DO_NOT_TOUCH);
  Matrix<T> SolveSystemUsingLU(Matrix<T> A,
                               Matrix<T> B,
                               const SolverMethod& method_ = SolverMethod::DO_NOT_TOUCH);
  Matrix<T> SolveSystem(Matrix<T> A,
                        Matrix<T> B,
                        const SolverMethod& method_ = SolverMethod::DO_NOT_TOUCH);
  Matrix<T> GetReversed(Matrix<T> A);
  Matrix<T> GetReversedAndDebugUsingDLUP(Matrix<T> A);

  void Print(const Matrix<T>& A, const vector<int>& ans_order) const;
  void Print(const Matrix<T>& A) const;
  void PrintThreeDiagonal(const vector<ThreeDiagonal<T>>& A) const;

  /// L D U P, such as A = D^{-1} L U P^{-1}
  tuple<Matrix<T>, Matrix<T>, Matrix<T>, Matrix<T>> DLUP_Step(Matrix<T>& A, int stage);
  tuple<Matrix<T>, Matrix<T>, Matrix<T>, Matrix<T>> DLUP_Decomposition(Matrix<T> A);

  tuple<Matrix<T>, vector<T>> LDL_Decomposition(Matrix<T> A);

  vector<T> SolveThreeDiagonalSystem(vector<ThreeDiagonal<T>> A, vector<T> b);

 private:

  void SolveStage(Matrix<T>& A, int stage, const SolverMethod& method_, vector<int>& ans_order);

  Matrix<T> makeExtended(const Matrix<T>& A, const Matrix<T>& B);
  Matrix<T> extractAnsMatrixFromExtended(const Matrix<T>& A);

  int best_in_the_row(const Matrix<T>& A, int row);
  int best_in_the_col(const Matrix<T>& A, int col);
  pair<int, int> best_in_the_sqr(const Matrix<T>& A, int start_i, int start_j);

  void swap_rows(Matrix<T>& A, int row1, int row2);
  void swap_columns(Matrix<T>& A, int col1, int col2, vector<int>& ans_order);
  void swap_columns(Matrix<T>& A, int col1, int col2);

  void normalize_row(Matrix<T>& A, int row);
  void substract_str(Matrix<T>& A, int row, int stage);
  void sub_row(Matrix<T>& A, int row1, int row2, T koef); // row1 -= row2 * koef
};


#endif //SYSTEM_SOLVER__SOLVER_H_
