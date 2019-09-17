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
      cout << j << " ";
    }
    cout << endl;
  }
  cout << "\n\n";
}


#endif //SYSTEM_SOLVER__DEBUG_H_
