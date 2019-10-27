//
// Created by user on 27.09.19.
//

#ifndef SYSTEM_SOLVER__HWSOLVER_H_
#define SYSTEM_SOLVER__HWSOLVER_H_

template<typename T>
class HW_Solver : public Solver<T> {
 public:
  Matrix<T> task1_crazy_gauss(Matrix<T> A) {
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
      // Solver<T>::sub_row(B, i - 1, i, k);
      for (int j = 0; j < n; ++j) {
        B[i - 1][j] -= B[i][j] * k;
      }
    }
    // PrintMatrix(A, "A");
    // PrintMatrix(B, "B");
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
        // Solver<T>::sub_row(B, i, j, A[i][j]);
        for (int k = 0; k < n; ++k) {
          B[i][k] -= B[j][k] * A[i][i];
          B[i][k] /= A[i][i];
        }
      }
      if (compareDouble(0, A[i][i])) throw logic_error("divide by zero");
      // A[i][i] = 1;
      // PrintMatrix(A, "A stage " + to_string(i));
      // PrintMatrix(B, "B stage " + to_string(i));
    }
    return B;
  }

  void task1_gauss_sub_row(Matrix<T>& A, Matrix<T>& B, int i, int n, int start_j, int end_j) {
    for (int j = start_j; j < end_j; ++j) {
      if (A[i][i] == 0) throw logic_error("divide by zero");
      T koef = A[j][i] / A[i][i];
      int mx = min(n, i + 2);
      for (int k = 0; k < i; ++k) {
        B[j][k] -= B[i][k] * koef;
      }
      for (int k = i; k < mx; ++k) {
        A[j][k] -= A[i][k] * koef;
        B[j][k] -= B[i][k] * koef;
      }
    }
  }

  Matrix<T> task1_gauss(Matrix<T>& A) {
    int n = A.size();
    Matrix < T > B(getIdentityMatrix<T>(n));
    for (int i = 0; i < n; ++i) {

      int left = i + 1, right = n;
      int mid = (left + right) / 2;

      std::future<void> handleFirst = async(launch::async, &HW_Solver::task1_gauss_sub_row,
                                            this, ref(A), ref(B), i, n, left, mid);
      std::future<void> handleSecond = async(launch::async, &HW_Solver::task1_gauss_sub_row,
                                             this, ref(A), ref(B), i, n, mid, right);
      //task1_gauss_sub_row(A, B, i, n, mid, right);
      handleFirst.get();
      handleSecond.get();
    }

    // PrintMatrix(A, "A end");

    for (int k = 0; k < n; ++k) {
      B.back()[k] /= A.back().back();
    }

    // PrintMatrix(B, "B end");

    for (int i = n - 2; i >= 0; --i) {
      for (int k = 0; k < n; ++k) {
        B[i][k] -= B[i + 1][k] * A[i][i + 1];
        B[i][k] /= A[i][i];
      }
    }
    return B;
  }

  Matrix<T> task1_random_strange_matrix(int n) {
    std::uniform_int_distribution<std::mt19937_64::result_type> udist(10, 1000);
    std::mt19937_64 generator(std::random_device{}());
    Matrix < T > ans(n, vector<T>(n));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < min(i + 2, n); ++j) {
        ans[i][j] = udist(generator);
      }
    }
    return ans;
  }

  T task5_get_Bij(int i, int j, int n) {
    if (!(i == 0 && j == 0) && !(i == n - 1 && j == n - 1) && (i == 0 || j == 0 || i == n - 1 || j == n - 1)) {
      return -1. / n;
    } else {
      return 0;
    }
  }


  /// ans and number of iterations

  pair<vector<double>, int> task5_relax(int n, double w = 1) {
    double EPS = 1e-10;
    vector<double> x(n, T(0));
    vector<double> diff(n), x_new;
    double error = 1;
    int iteration_count = 0;
    int iteration_limit = 100;
    while (error > EPS && iteration_count < iteration_limit) {
      //cout << error << endl;
      ++iteration_count;
      x_new = task5_relax_step(x, w);
      for (int i = 0; i < n; ++i) {
        diff[i] = abs(x[i] - x_new[i]);
      }
      error = *max_element(begin(diff), end(diff));
      x = x_new;
    }
    return {x, iteration_count};
  }

  vector<double> task5_relax_step(const vector<double>& x, double w) {
    int n = x.size();
    vector<double> new_x(n, 1. / n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
        new_x[i] += task5_get_Bij(i, j, n) * new_x[j];
      }
      for (int j = i; j < n; ++j) {
        new_x[i] += task5_get_Bij(i, j, n) * x[j];
      }
      new_x[i] *= w;
      new_x[i] += (1. - w) * x[i];
    }
    return new_x;
  }

};

#endif //SYSTEM_SOLVER__HWSOLVER_H_
