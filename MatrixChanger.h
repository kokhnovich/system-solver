//
// Created by user on 26.10.19.
//

#ifndef SYSTEM_SOLVER__MATRIXCHANGER_H_
#define SYSTEM_SOLVER__MATRIXCHANGER_H_

// template<class T>
class MatrixChanger {
 public:
  MatrixChanger(const Matrix<Fraction>& a) : a_(a) {}

  /*
   * SwapRows initializes SwapColumns and vice versa
  */

  void SwapRows(int i, int j, bool and_cols = true) {
    a_[i].swap(a_[j]);
    cerr << "swap rows: " + to_string(i) + " " + to_string(j) << endl;
    if (and_cols) {
      SwapColumns(j, i, false);
    }
  }

  void SwapColumns(int i, int j, bool and_rows = true) {
    for (int k = 0; k < a_.size(); ++k) {
      swap(a_[k][i], a_[k][j]);
    }
    cerr << "swap columns: " + to_string(i) + " " + to_string(j) << endl;
    if (and_rows) {
      SwapRows(j, i, false);
    }
  }

  void MultRow(int i, Fraction koef, bool first = true) {
    for (int j = 0; j < a_[i].size(); ++j) {
      a_[i][j] *= koef;
    }
    cerr << "mult row: " << i << " " << koef << endl;
    if (first) {
      MultCol(i, Fraction(1) / koef, false);
    }
  }

  void MultCol(int j, Fraction koef, bool first = true) {
    for (int i = 0; i < a_.size(); ++i) {
      a_[i][j] *= koef;
    }
    cerr << "mult col: " << j << " " << koef << endl;
    if (first) {
      MultRow(j, Fraction(1) / koef, false);
    }
  }

  void AddRow(int i, int j, Fraction koef, bool first = true) {
    for (int k = 0; k < a_[i].size(); ++k) {
      a_[i][k] += a_[j][k] * koef;
    }
    cerr << "add row: " << i << " " << j << " " << koef << endl;
    if (first) {
      AddCol(j, i, koef, false);
    }
  }

  void AddCol(int i, int j, Fraction koef, bool first = true) {
    for (auto& k : a_) {
      k[i] += k[j] * koef;
    }
    cerr << "add col: " << i << " " << j << " " << koef << endl;
    if (first) {
      AddRow(j, i, koef, false);
    }
  }

  void Print() {
    PrintMatrix(a_, "A");
  }

 private:
  Matrix<Fraction> a_;
};

#endif //SYSTEM_SOLVER__MATRIXCHANGER_H_
