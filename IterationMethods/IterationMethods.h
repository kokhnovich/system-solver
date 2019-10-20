
template<typename T>
class IterationMethods {
 public:

  vector<T> YakobyMethod(Matrix<T> A, vector<T> b) {
    int n = A.size();
    Matrix < T > B(n, vector<T>(n, T(0)));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i != j) {
          B[i][j] = -A[i][j] / A[i][i];
        }
      }
    }

    PrintMatrix(B, "B");

    // @TODO make assertion

    vector<T> d(n);
    for (int i = 0; i < n; ++i) {
      d[i] = b[i] / A[i][i];
    }

    vector<T> prev_x(n, T(0));

    int iters = 10;
    while (iters--) {
      PrintMatrix(prev_x, "iter " + to_string(iters));
      prev_x = YakobyMethodStep(prev_x, B, d);
    }
    return prev_x;
  }

  vector<T> GaussZeydelMethod(Matrix<T> A, vector<T> b) {
    int n = A.size();
    Matrix < T > B(n, vector<T>(n, T(0)));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i != j) {
          B[i][j] = -A[i][j] / A[i][i];
        }
      }
    }
    PrintMatrix(B, "B");
    // @TODO make assertion
    vector<T> d(n);
    for (int i = 0; i < n; ++i) {
      d[i] = b[i] / A[i][i];
    }
    vector<T> prev_x(n, T(0));
    int iters = 10;
    while (iters--) {
      PrintMatrix(prev_x, "iter " + to_string(iters));
      prev_x = GaussZeydelMethodStep(prev_x, B, d);
    }
    return prev_x;
  }

  
  vector<T> RelaxationMethod(Matrix<T> A, vector<T> b, double w = 0.5) {
    int n = A.size();
    Matrix < T > B(n, vector<T>(n, T(0)));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i != j) {
          B[i][j] = -A[i][j] / A[i][i];
        }
      }
    }
    PrintMatrix(B, "B");
    // @TODO make assertion
    vector<T> d(n);
    for (int i = 0; i < n; ++i) {
      d[i] = b[i] / A[i][i];
    }
    vector<T> prev_x(n, T(0));
    int iters = 10;
    while (iters--) {
      PrintMatrix(prev_x, "iter " + to_string(iters));
      prev_x = RelaxationMethodStep(prev_x, B, d, w);
    }
    return prev_x;
  }

 private:
  vector<T> YakobyMethodStep(const vector<T>& x, const Matrix<T>& B, const vector<T>& d) {
    int n = B.size();
    vector<T> new_x(d);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        new_x[i] += B[i][j] * x[j];
      }
    }
    return new_x;
  }
  vector<T> GaussZeydelMethodStep(const vector<T>& x, const Matrix<T>& B, const vector<T>& d) {
    int n = B.size();
    vector<T> new_x(d);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
        new_x[i] += B[i][j] * new_x[j];
      }
      for (int j = i; j < n; ++j) {
        new_x[i] += B[i][j] * x[j];
      }
    }
    return new_x;
  }
  vector<T> RelaxationMethodStep(const vector<T>& x, const Matrix<T>& B, const vector<T>& d, double w) {
    int n = B.size();
    vector<T> new_x(d);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
        new_x[i] += B[i][j] * new_x[j];
      }
      for (int j = i; j < n; ++j) {
        new_x[i] += B[i][j] * x[j];
      }
      new_x[i] *= w;
      new_x[i] += (1. - w) * x[i];
    }
    return new_x;
  }

};