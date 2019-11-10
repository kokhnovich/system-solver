//
// Created by user on 27.09.19.
//

#ifndef SYSTEM_SOLVER__HWSOLVER_H_
#define SYSTEM_SOLVER__HWSOLVER_H_

template<typename T>
class HW_Solver : public Solver<T> {
 public:
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

  Matrix<T> task1_gauss(Matrix<T> A) {
    int n = A.size();
    Matrix<T> B(getIdentityMatrix<T>(n));
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

  // don't store zeros
  Matrix<T> task1_random_strange_matrix(int n) {
    std::uniform_int_distribution<std::mt19937_64::result_type> udist(10, 1000);
    std::mt19937_64 generator(std::random_device{}());
    Matrix<T> ans(n, vector<T>(n));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < min(i + 2, n); ++j) {
        ans[i][j] = udist(generator);
      }
    }
    return ans;
  }

  vector<double> task2_dlup(Matrix<double>& U, vector<double>& b) {
    int n = U.size();
    Matrix<double> L(n, vector<T>(n, 0));
    vector<pair<int, int>> d;
    bool is_first = true;
    // alongside finding L^T we can do the same operations with vector b
    // instead of solving lower triangular matrix
    for (int stage = 0; stage < n; ++stage) {
      pair<int, int> best_pos = {Solver<double>::best_in_the_col(U, stage), stage};

      Solver<double>::swap_rows(U, stage, best_pos.first);
      Solver<double>::swap_columns(U, stage, best_pos.second);

      d.push_back(std::make_pair(stage, best_pos.first));

      if (is_first) {
        is_first = false;
      } else {
        Solver<double>::swap_rows(L, stage, best_pos.first);
      }

      for (int row = stage; row < n; ++row) {
        L[row][stage] = U[row][stage];
      }

      assert(U[stage][stage] != T(0));
      double k = U[stage][stage];
      for (int col = stage; col < n; ++col) {
        U[stage][col] /= k;
      }

      for (int row = stage + 1; row < n; ++row) {
        Solver<double>::sub_row(U, row, stage, U[row][stage]);
      }
    }

    vector<int> d_order(n);
    for (int i = 0; i < n; ++i) {
      d_order[i] = i;
    }

    for (int i = d.size() - 1; i >= 0; --i) {
      std::swap(d_order[d[i].first], d_order[d[i].second]);
    }

    vector<double> new_b(n);
    for (int i = 0; i < n; ++i) {
      new_b[d_order[i]] = b[i];
    }
    b = (new_b);

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

  vector<double> task3_lu_for_sym(Matrix<double>& U, vector<double>& b) {
    int n = U.size();
    for (int i = 1; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
        U[i][j] = 0;
      }
    }

    // do the same operations with vector b instead of solving a system with lower triangular matrix
    Matrix<double> L(n, vector<double>(n, 0));
    for (int stage = 0; stage < n; ++stage) {
      L[stage][stage] = 1;
      for (int row = stage + 1; row < n; ++row) {
        double k = U[stage][row] / U[stage][stage];
        L[row][stage] = k;
        for (int col = row; col < n; ++col) {
          U[row][col] -= k * U[stage][col];
        }
        b[row] -= b[stage] * k;
      }
    }

    return Solver<double>::SolveUp(U, b);
  }

  vector<double> task3_ldlt(Matrix<double>& L, vector<double>& b) {
    // vector<T> D(L.size(), 1);
    // instead of storing D I just multiply elements from vector b by -1
    int n = L.size();

    for (int now = 0; now < n; ++now) {
      double k = (L[now][now] < 0 ? -sqrt(-L[now][now]) : sqrt(L[now][now]));
      for (int col = now; col < n; ++col) {
        L[now][col] /= k;
      }
      b[now] /= k;
      for (int row = now + 1; row < n; ++row) {
        double kk = L[row][now] / L[now][now];
        for (int col = now; col < L.size(); ++col) {
          L[row][col] -= L[now][col] * kk;
        }
        b[row] -= kk * b[now];
      }
    }

    // instead of transposition we can use L[j][i] rather than L[i][j] in the code below
    // L = Solver<double>::Transpose(L);

    // PrintMatrix(mult(L, Solver<double>::Transpose(L)));
    // PrintMatrix(b, "b 2.0");

    for (int i = n - 1; i >= 0; --i) {
      for (int j = i + 1; j < n; ++j) {
        b[i] -= L[i][j] * b[j];
      }
      b[i] /= L[i][i];
    }
    return b;
  }

  vector<double> task4_threediagonal(Matrix<double> A, vector<double> b) {
    int n = A.size();
    vector<double> d(n, 0);
    for (int i = 1; i < n; ++i) {
      if (A[i - 1][1] == 0) throw std::logic_error("divide by zero");
      double k = A[i][0] / A[i - 1][1];

      A[i][0] -= k * A[i - 1][1];
      A[i][1] -= k * A[i - 1][2];
      b[i] -= k * b[i - 1];
    }

    // PrintThreeDiagonal(A);
    // PrintMatrix(b, "b");

    b.back() /= A.back()[1];
    A.back()[1] = 1;
    for (int i = n - 2; i >= 0; --i) {
      b[i] -= A[i][2] * b[i + 1];
      b[i] /= A[i][1];
      A[i][1] = 1;
    }
    return b;
  }

  T task5_get_Bij(int i, int j, int n) {
    if (!(i == 0 && j == 0) && !(i == n - 1 && j == n - 1) && (i == 0 || j == 0 || i == n - 1 || j == n - 1)) {
      return -1. / n;
    } else {
      return 0;
    }
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
  /// ans and number of iterations
  pair<vector<double>, int> task5_relax(int n, double w = 1) {
    double EPS = 1e-10;
    // don't store matrix at all: function task5_get_Bij returns B[i][j] in O(const)
    vector<double> x(n, T(0));
    vector<double> diff(n), x_new;
    double error = 1;
    int iteration_count = 0;
    int iteration_limit = 100; // don't wait if it tooo slooow
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
};

#endif //SYSTEM_SOLVER__HWSOLVER_H_
