template<typename T>
Matrix<T> mult(const Matrix<T>& A, const Matrix<T>& B) {
  if (A.empty() || B.empty() || A[0].size() != B.size()) throw logic_error("wrong sizes");
  size_t N = A.size(), M = A[0].size(), K = B[0].size();
  Matrix<T> R(N, vector<T>(K, T(0)));
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < K; j++) {
      for (size_t k = 0; k < M; k++) {
        R[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return R;
}

template<typename T>
Matrix<T> mult(const Matrix<T>& A, const vector<T>& B) {
  return mult(A, {B});
}

template<typename T>
Matrix<T> mult(const vector<T>& A, const Matrix<T>& B) {
  return mult({A}, B);
}

template<typename T>
Matrix<T> mult(const vector<T>& A, const vector<T>& B) {
  return mult({A}, {B});
}

template<typename T>
Matrix<T> Transpose(Matrix<T> A) {
  assert(A.size() == A[0].size());
  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      swap(A[i][j], A[j][i]);
    }
  }
  return A;
}

void GetMatrixFromSTDIN(Matrix<double>& a) {
  for (int i = 0; i < a.size(); ++i) {
    for (int j = 0; j < a[i].size(); ++j) {
      cin >> a[i][j];
    }
  }
}

template<typename T>
Matrix<T> generateThreeDiagonal(int size_, int min_elem = -100, int max_elem = 100) {
  std::uniform_int_distribution<std::mt19937_64::result_type> udist(min_elem, max_elem);
  std::mt19937_64 generator(std::random_device{}());
  Matrix<T> a(size_, vector<T>(3));
  for (auto& i : a) {
    for (auto& j : i) {
      j = udist(generator);
    }
  }
  return a;
}

template<typename T>
Matrix<T> generateRandomSymmetricMatrix(int size_, int min_elem = -100, int max_elem = 100) {
  std::uniform_int_distribution<std::mt19937_64::result_type> udist(min_elem, max_elem);
  std::mt19937_64 generator(std::random_device{}());
  Matrix<T> a(size_, vector<T>(size_));
  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < size_; ++j) {
      if (j < i) {
        a[i][j] = a[j][i];
      } else {
        a[i][j] = udist(generator);
      }
    }
  }
  return a;
}

double EPS = 1e-3;

inline bool compareDouble(double d1, double d2) {
  return abs(d1 - d2) <= EPS;
}

inline bool compareVectorOfDouble(const vector<double>& A1, const vector<double>& A2) {
  assert(A1.size() == A2.size());
  for (int i = 0; i < A1.size(); ++i) {
    if (!compareDouble(A1[i], A2[i])) {
      return false;
    }
  }
  return true;
}

inline bool compareMatrixOfDouble(const Matrix<double>& A1, const Matrix<double>& A2) {
  assert(A1.size() == A2.size());
  for (int i = 0; i < A1.size(); ++i) {
    if (!compareVectorOfDouble(A1[i], A2[i])) {
      return false;
    }
  }
  return true;
}

template<typename T>
vector<T> generateAnsMatrix(const Matrix<T>& A, const vector<T>& x) {
  vector<T> ans(A.size());
  for (int i = 0; i < A.size(); ++i) {
    ans[i] = T(0);
    for (int j = 0; j < A.size(); ++j) {
      ans[i] += x[j] * A[i][j];
    }
  }
  return ans;
}

/// Nx1 == true  <=> ans.size() == n && ans[i].size() == 1 \forall i
/// Nx1 == false <=> ans.size() == 1 && ans[i].size() == n \forall i
template<typename T>
Matrix<T> make2Dfrom1D(const vector<T>& a, bool Nx1 = true) {
  if (Nx1) {
    Matrix<T> ans(a.size(), vector<T>(1, 0));
    for (int i = 0; i < a.size(); ++i) {
      ans[i][0] = a[i];
    }
    return ans;
  } else {
    Matrix<T> ans(1, vector<T>(a.size()));
    for (int i = 0; i < a.size(); ++i) {
      ans[0][i] = a[i];
    }
    return ans;
  }
}

template<typename T>
vector<T> make1Dfrom2D(const Matrix<T>& a) {
  if (a.size() != 1) {
    vector<T> ans(a.size());
    for (int i = 0; i < a.size(); ++i) {
      ans[i] = a[i][0];
    }
    return ans;
  } else {
    vector<T> ans(a[0].size());
    for (int i = 0; i < a[0].size(); ++i) {
      ans[i] = a[0][i];
    }
    return ans;
  }
}

template<typename T>
vector<T> make1Dfrom2DForOrderMatrixes(const Matrix<T>& a) {
  vector<T> ans(a.size());
  for (int i = 0; i < a.size(); ++i) {
    bool ok = false;
    for (int j = 0; j < a.size(); ++j) {
      if (a[i][j]) {
        ans[i] = j;
        ok = true;
        break;
      }
    }
    if (!ok) {
      assert(false);
    }
  }
  return ans;
}

template<class T>
struct greater_using_abs {
  bool operator()(const T& x, const T& y) const { return abs(x) > abs(y); }
};

template<typename T>
Matrix<T> getHWMatrix(int n) {
  return {{n + 1, n / 2, -n / 2, 1},
          {-n - 1, -(n + 1) / 2, (n + 1) / 2, -(n + 2) / 3},
          {-n + 1, (n + 1) / 2, -(n + 2) / 3, n - 1},
          {n / 3, -1, n, -n}};
}

Matrix<double> getNewHWMatrix(double n = 12) {
  Matrix<double> m = {{2 - n, 0, n - 1},
                      {n - 1, 1, 1 - n},
                      {-2 * (n - 1), 0, 2 * n - 1}};
  return m;
}

Matrix<double> getAnotherHWMatrix(double n = 12) {
  Matrix<double> m = {{2 * (n + 1), 4 + 5 * n, 1 + 2 * n},
                      {-1 - 2 * n, -3 - 5 * n, -1 - 2 * n},
                      {3 * n, 7 * n + 2, 1 + 3 * n}};
  return m;
}

template<typename T>
Matrix<T> getIdentityMatrix(int n) {
  Matrix<T> A(n, vector<T>(n, T(0)));
  for (size_t i = 0; i < n; ++i) {
    A[i][i] = T(1);
  }
  return A;
}

template<typename T>
T scalar_mult_vectors(const vector<T>& a, const vector<T>& b) {
  T ans = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    ans += a[i] * b[i];
  }
}

template<typename T>
vector<T> operator+(vector<T> a, const vector<T>& b) {
  if (a.size() != b.size()) {
    throw logic_error("vectors are different to sub");
  } else {
    for (int i = 0; i < a.size(); ++i) {
      a[i] += b[i];
    }
    return a;
  }
}

template<typename T>
vector<T> operator-(vector<T> a, const vector<T>& b) {
  if (a.size() != b.size()) {
    throw logic_error("vectors are different to sub");
  } else {
    for (int i = 0; i < a.size(); ++i) {
      a[i] -= b[i];
    }
    return a;
  }
}

double GetMatrixNorm(const vector<double>& a) {
  double res = 0;
  for (auto& i : a) {
    res += i;
  }
  return (double) (res / a.size());
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b) {
  return mult(a, b);
}

template<typename T>
ostream& operator<<(ostream& os, const Matrix<T>& a) {
  for(auto& i : a) {
    for(auto& j : i) {
      os << j << " ";
    }
    os << endl;
  }
  os << endl;
}

template<typename T>
ostream& operator<<(ostream& os, const vector<T>& a) {
  for(auto& i : a) {
    os << i << " ";
  }
  os << endl;
}