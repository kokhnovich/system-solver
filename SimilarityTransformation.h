//
// Created by user on 1.12.19.
//

#ifndef SYSTEM_SOLVER__SIMILARITYTRANSFORMATION_H_
#define SYSTEM_SOLVER__SIMILARITYTRANSFORMATION_H_

// #define SimTransDebug

namespace SimilarityTransformation {

void SwapRows(Matrix<double>& A, int i, int j, bool and_cols = true);
void SwapCols(Matrix<double>& A, int i, int j, bool and_cols = true);
void MultRow(Matrix<double>& A, int i, double koef, bool first = true);
void MultCol(Matrix<double>& A, int i, double koef, bool first = true);
void AddRow(Matrix<double>& A, int i, int j, double koef, bool first = true);
void AddCol(Matrix<double>& A, int i, int j, double koef, bool first = true);

void SwapRows(Matrix<double>& A, int i, int j, bool and_cols) {
  A[i].swap(A[j]);
  //cerr << "swap rows: " + to_string(i) + " " + to_string(j) << endl;
  if (and_cols) {
    SwapCols(A, j, i, false);
  }
}

void SwapCols(Matrix<double>& A, int i, int j, bool and_rows) {
  for (auto& k : A) {
    swap(k[i], k[j]);
  }
  //cerr << "swap columns: " + to_string(i) + " " + to_string(j) << endl;
  if (and_rows) {
    SwapRows(A, j, i, false);
  }
}

void MultRow(Matrix<double>& A, int i, double koef, bool first) {
  for (int j = 0; j < A[i].size(); ++j) {
    A[i][j] *= koef;
  }
  //cerr << "mult row: " << i << " " << koef << endl;
  if (first) {
    MultCol(A, i, double(1) / koef, false);
  }
}

void MultCol(Matrix<double>& A, int j, double koef, bool first) {
  for (int i = 0; i < A.size(); ++i) {
    A[i][j] *= koef;
  }
  //cerr << "mult col: " << j << " " << koef << endl;
  if (first) {
    MultRow(A, j, double(1) / koef, false);
  }
}

void AddRow(Matrix<double>& A, int i, int j, double koef, bool first) {
  for (int k = 0; k < A[i].size(); ++k) {
    A[i][k] += A[j][k] * koef;
  }
  //cerr << "add row: " << i << " " << j << " " << koef << endl;
  if (first) {
     AddCol(A, j, i, -koef, false);
  }
}

void AddCol(Matrix<double>& A, int i, int j, double koef, bool first) {
  for (auto& k : A) {
    k[i] += k[j] * koef;
  }
  //cerr << "add col: " << i << " " << j << " " << koef << endl;
  if (first) {
    AddRow(A, j, i, -koef, false);
  }
}


}

#endif //SYSTEM_SOLVER__SIMILARITYTRANSFORMATION_H_
