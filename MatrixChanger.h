//
// Created by user on 26.10.19.
//

#ifndef SYSTEM_SOLVER__MATRIXCHANGER_H_
#define SYSTEM_SOLVER__MATRIXCHANGER_H_

// template<class T>
class MatrixChanger {
 public:
  MatrixChanger(const Matrix<double>& a) : a_(a) {}

  /*
   * SwapRows initializes SwapColumns and vice versa
  */

  void SwapRows(int i, int j, bool and_cols = true) {
    a_[i].swap(a_[j]);
    cerr << "swap rows: " + to_string(i) + " " + to_string(j) << endl;
    if (and_cols) {
      SwapCols(j, i, false);
    }
  }

  void SwapCols(int i, int j, bool and_rows = true) {
    for (auto& k : a_) {
      swap(k[i], k[j]);
    }
    cerr << "swap columns: " + to_string(i) + " " + to_string(j) << endl;
    if (and_rows) {
      SwapRows(j, i, false);
    }
  }

  void MultRow(int i, double koef, bool first = true) {
    for (int j = 0; j < a_[i].size(); ++j) {
      a_[i][j] *= koef;
    }
    cerr << "mult row: " << i << " " << koef << endl;
    if (first) {
      MultCol(i, double(1) / koef, false);
    }
  }

  void MultCol(int j, double koef, bool first = true) {
    for (int i = 0; i < a_.size(); ++i) {
      a_[i][j] *= koef;
    }
    cerr << "mult col: " << j << " " << koef << endl;
    if (first) {
      MultRow(j, double(1) / koef, false);
    }
  }

  void AddRow(int i, int j, double koef, bool first = true) {
    for (int k = 0; k < a_[i].size(); ++k) {
      a_[i][k] += a_[j][k] * koef;
    }
    cerr << "add row: " << i << " " << j << " " << koef << endl;
    if (first) {
      AddCol(j, i, -koef, false);
    }
  }

  void AddCol(int i, int j, double koef, bool first = true) {
    for (auto& k : a_) {
      k[i] += k[j] * koef;
    }
    cerr << "add col: " << i << " " << j << " " << koef << endl;
    if (first) {
      AddRow(j, i, -koef, false);
    }
  }

  void Print() {
    PrintMatrix(a_, "A");
  }

 private:
  Matrix<double> a_;
};

#endif //SYSTEM_SOLVER__MATRIXCHANGER_H_
