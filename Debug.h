//
// Created by user on 13.09.19.
//

#ifndef SYSTEM_SOLVER__DEBUG_H_
#define SYSTEM_SOLVER__DEBUG_H_

template<typename T>
void PrintMatrix(const vector<vector<T>>& a, const string& message = "") {
  cout << message << ":\n";
  for (auto& i : a) {
    for (auto& j : i) {
      cout << setw(6) << j << " ";
    }
    cout << endl;
  }
  cout << "\n\n";
}
template<typename T>
void PrintMatrix(const vector<T>& a, const string& message = "") {
  cout << message << ":\n";
  for (auto& i : a) {
    cout << i << " ";
  }
  cout << "\n\n";
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
  cout << "\n\n";
}

#endif //SYSTEM_SOLVER__DEBUG_H_
