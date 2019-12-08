//
// Created by user on 25.11.19.
//

#ifndef SYSTEM_SOLVER_ITERATIONMETHODS_EIGENPAIRSFINDITERATIONMETHODS_H_
#define SYSTEM_SOLVER_ITERATIONMETHODS_EIGENPAIRSFINDITERATIONMETHODS_H_

double sq(double a) {
  return a * a;
}

double getNorm(const vector<double>& a) {
  double ans = -1;
  for (const auto& i : a) {
    ans = max(ans, abs(i));
  }
  return ans;
}

// vector of eigenpairs - double and vector<double>
vector<pair<double, vector<double>>> getMaximalEigenValuesWithTheirEigenvectors(Matrix<double> A) {
  int n = A.size();

  vector<double> prev(n, 0.);
//  for (int i = 0; i < n; i += 2) {
//    prev[i] = 1.;
//  }
  prev[0] = 1.;

  double norm_prev = n / 2;
  double res_prev = 0;

  for (int iter = 0; iter < 1000; ++iter) {
    vector<double> now(n, 0.);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        now[i] += A[i][j] * prev[j];
      }
    }

    double norm = getNorm(now);

    for (auto& i : now) {
      i /= norm;
    }

//    if (iter % 1 == 0) {
    double res = sqrt(norm * norm_prev);
    cout << norm << " " << res << " " << abs(res - res_prev) << endl;
    res_prev = res;
//    }
    PrintMatrix(now, to_string(iter));
    norm_prev = norm;
    prev = now;
  }

  return {};
};

vector<double> applyDanilevskyMethod(Matrix<double> A) {
  int n = A.size();
  vector<double> polynom;
  polynom.reserve(n + 1);
  for (int stage = n - 2; stage >= 0; --stage) {

    int col = stage;
    int row = stage + 1;

    bool found = false;
    for (int j = col; j >= 0; --j) {
      if (A[j + 1][j] != 0) {
        SimilarityTransformation::SwapCols(A, j, col);
        found = true;
        break;
      }
    }

    if (!found) throw std::invalid_argument("bad matrix");

    auto M = getIdentityMatrix<double>(n);
    for (int j = 0; j < n; ++j) {
      M[col][j] = -A[row][j] / A[row][col];
    }
    M[col][col] = 1. / A[row][col];

    auto M_rev = getIdentityMatrix<double>(n);
    for (int j = 0; j < n; ++j) {
      M_rev[col][j] = A[row][j];
    }

    //PrintMatrix(M, "M");
    //PrintMatrix(M_rev, "M_rev");
    A = M_rev * A * M;
    //PrintMatrix(A, "M_rev A M");
  }
  polynom.push_back(1.);
  for (auto& i : A[0]) {
    if (abs(i) > EPS) {
      polynom.push_back(-i);
    } else {
      polynom.push_back(0.);
    }
  }
  return polynom;
}

vector<double> findRootsFromRangeBruteForce(const vector<double>& a, double mn, double mx, double step) {
  vector<double> roots;
  for (double root = mn; root <= mx; root += step) {
    double value = 0;
    double lambda = 1;
    for (int i = a.size() - 1; i >= 0; --i) {
      value += lambda * a[i];
      lambda *= root;
    }
    if (abs(value) < 1.) {
      roots.push_back(root);
    }
  }
  return roots;
}

void QRDecomposition(Matrix<double> A) {

}

void findingEigenpairsUsingQRAlgorithmWithReflectionsMethod(Matrix<double> A) {
  int n = A.size();

  Matrix<double> Q(getIdentityMatrix<double>(n));
  Matrix<double> A_copy(A);

  PrintMatrix(A, "A start");
  for (int j = 0; j < n - 1; ++j) {
    for (int i = n - 1; i > j; --i) {
      double cos = A[j][j] / (sq(A[i][j]) + sq(A[j][j]));
      double sin = A[i][j] / (sq(A[i][j]) + sq(A[j][j]));

      vector<double> new_i = A[i];
      vector<double> new_j = A[j];

      for (int col = 0; col < n; ++col) {
        new_i[col] = cos * A[i][col] - sin * A[j][col];
        new_j[col] = sin * A[i][col] + cos * A[j][col];
      }

      A[i].swap(new_i);
      A[j].swap(new_j);

      Matrix<double> Q_local(getIdentityMatrix<double>(n));
      Q_local[j][j] = cos;
      Q_local[j][i] = sin;
      Q_local[i][j] = -sin;
      Q_local[i][i] = cos;

      Q = Q_local * Q;

      PrintMatrix(Q_local, "Q_local");
      PrintMatrix(A, "A");
      PrintMatrix(Q_local*A_copy, "QA");
    }
  }
  PrintMatrix(Q, "Q");
  PrintMatrix(Q*A_copy, "QR");

  // A
  // A_copy = Q^T * A

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if (abs(A[i][j]) < EPS) {
        A[i][j] = 0.;
      }
    }
  }
  PrintMatrix(A);

}

#endif //SYSTEM_SOLVER_ITERATIONMETHODS_EIGENPAIRSFINDITERATIONMETHODS_H_
