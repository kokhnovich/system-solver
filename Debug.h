//
// Created by user on 13.09.19.
//

#ifndef SYSTEM_SOLVER__DEBUG_H_
#define SYSTEM_SOLVER__DEBUG_H_

template<typename T>
void PrintMatrixInLatex(const vector<vector<T>>& a, const string& message = "") {
  cout << message << endl;
  cout << "\\begin{bmatrix}" << endl;
  for (auto& i : a) {
    for (auto& j : i) {
      cout << setw(6) << j << " & ";
    }
    cout << "\\\\" << endl;
  }
  cout << "\\end{bmatrix}" << endl;
}

template<typename T>
void PrintMatrix(const vector<vector<T>>& a, const string& message = "", bool latex=LATEX) {
  if (latex) {
    PrintMatrixInLatex(a, message);
  } else {
    cout << message << ":\n";
    for (auto& i : a) {
      for (auto& j : i) {
        cout << setw(6) << j << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}

template<typename T>
void PrintMatrixInLatex(const vector<T>& a, const string& message = "") {
  cout << message << endl;
  cout << "\\begin{bmatrix}" << endl;

  for (auto& j : a) {
    cout << setw(6) << j << " & ";
  }
  cout << "\\\\" << endl;

  cout << "\\end{bmatrix}" << endl;
}

template<typename T>
void PrintMatrix(const vector<T>& a, const string& message = "", bool latex=LATEX) {
  if (latex) {
    PrintMatrixInLatex(a, message);
  } else {
    cout << message << ":\n";
    for (auto& i : a) {
      cout << i << " ";
    }
    cout << endl;
  }
}


template<typename T>
void PrintPartOfMatrix(const vector<vector<T>>& a, size_t size_ = 8, const string& message = "") {
  cout << message << ":\n";
  for (int i = 0; i < min(size_, a.size()); ++i) {
    for (int j = 0; j < min(size_, a[i].size()); ++j) {
      cout << a[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

#endif //SYSTEM_SOLVER__DEBUG_H_
