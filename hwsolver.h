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

  vector<double> task2_dlup(Matrix<double>& A, vector<double>& b) {
    int n = A.size();

    vector<Matrix<T>> Ps, Ds;
    Matrix<double> L(n, vector<T>(n, 0)), U(move(A)), I(getIdentityMatrix<T>(n));

    vector<pair<int, int>> p, d;

    bool is_first = true;
    for (int stage = 0; stage < n; ++stage) {
      pair<int, int> best_pos = {Solver<double>::best_in_the_col(U, stage), stage};

      Solver<double>::swap_rows(U, stage, best_pos.first);
      Solver<double>::swap_columns(U, stage, best_pos.second);

//    if (stage != best_pos.first) {
//      cout << "change rows fact " << stage << " " << best_pos.first << "\n";
//    }
//    if (stage != best_pos.second) {
//      cout << "change column fact " << stage << " " << best_pos.second << "\n";
//    }

      Matrix<double> P(I), D(I);
      swap(P[stage], P[best_pos.second]);
      swap(D[stage], D[best_pos.first]);
      Ds.push_back(D);
      Ps.push_back(P);

      p.push_back(std::make_pair(stage, best_pos.second));
      d.push_back(std::make_pair(stage, best_pos.first));

      // PrintMatrix(U, "U after changes");

      if (is_first) {
        is_first = false;
      } else {
        L = mult(D, L);
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

//    cout << "After stage " << stage << endl;
//    PrintMatrix(U, "U" + to_string(stage));
//    PrintMatrix(L, "L" + to_string(stage));
//    PrintMatrix(P, "P" + to_string(stage));
//    PrintMatrix(D, "D" + to_string(stage));
    }

    Matrix<double> P = I;
    for (int i = 0; i < Ps.size(); ++i) {
      P = mult(P, Ps[i]);
    }

    vector<int> p_order(n), d_order(n);
    for(int i = 0; i < n; ++i) {
      p_order[i] = i;
      d_order[i] = i;
    }

    Matrix<double> D = I;
    for (int i = Ds.size() - 1; i >= 0; --i) {
      D = mult(D, Ds[i]);
    }

    for(int i = d.size() - 1; i >= 0; --i) {
      std::swap(d_order[d[i].first], d_order[d[i].second]);
    }
    for(int i = 0; i < p.size(); ++i) {
      std::swap(p_order[p[i].first], p_order[p[i].second]);
    }

    PrintMatrix(p_order, "p_order");
    PrintMatrix(d_order, "d_order");

    PrintMatrix(D, "D");
    PrintMatrix(L, "L");
    PrintMatrix(U, "U");
    PrintMatrix(P, "P");

    vector<double> b_original(b);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (D[i][j]) {
          b[i] = b_original[j];
        }
      }
    }
    // L Y = B
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
