//
// Created by user on 27.09.19.
//

#ifndef SYSTEM_SOLVER__HWSOLVER_H_
#define SYSTEM_SOLVER__HWSOLVER_H_

// #include "Solver.cpp"

template<typename T>

class HW_Solver : public Solver<T> {
 public:
  Matrix<T> task2_gauss(Matrix<T> A) {
    int n = A.size();
    Matrix < T > B(getIdentityMatrix<T>(n));
    for (int i = n - 1; i >= 1; --i) {
      // cout << i << " " << A[i - 1][i] / A[i][i] << endl;
      if (compareDouble(A[i][i], 0)) throw logic_error("divide by zero");
      T k = A[i - 1][i] / A[i][i];
      // Solver<T>::sub_row(A, i - 1, i, k);
      for (int j = 0; j <= i; ++j) {
        A[i - 1][j] -= A[i][j] * k;
      }
      Solver<T>::sub_row(B, i - 1, i, k);
    }
    // PrintMatrix(A, "A");
    // PrintMatrix(B, "B");
    for (int i = 0; i < A.size(); ++i) {
      for (int j = 0; j < i; ++j) {
        Solver<T>::sub_row(B, i, j, A[i][j]);
      }
      if (compareDouble(0, A[i][i])) throw logic_error("divide by zero");
      for (int j = 0; j < B[i].size(); ++j) {
        B[i][j] /= A[i][i];
      }
      // A[i][i] = 1;
      // PrintMatrix(A, "A stage " + to_string(i));
      // PrintMatrix(B, "B stage " + to_string(i));
    }
    return B;
  }

  Matrix<T> task2_random_strange_matrix(int n) {
    std::uniform_int_distribution<std::mt19937_64::result_type> udist(10050, 1005000);
    std::mt19937_64 generator(std::random_device{}());
    Matrix < T > ans(n, vector<T>(n));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < min(i + 2, n); ++j) {
        ans[i][j] = udist(generator);
      }
    }
    return ans;
  }

};

#endif //SYSTEM_SOLVER__HWSOLVER_H_
