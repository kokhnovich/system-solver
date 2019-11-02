/*tuple<Matrix<double>, Matrix<double>, Matrix<double>, Matrix<double>> DLUP_Decomposition(Matrix<double> A,
                                                                                         const SolverMethod& method_) {
  int n = A.size();
  vector<Matrix<double>> Ps, Ds;
  Matrix<double> L(n, vector<double>(n, T(0))), U(move(A)), I(getIdentityMatrix<T>(n));

  bool is_first = true;
  for (int stage = 0; stage < n; ++stage) {
    pair<int, int> best_pos;
    switch (method_) {
      case SolverMethod::BEST_IN_MATRIX: {
        best_pos = best_in_the_sqr(U, stage, stage);
        break;
      }
      case SolverMethod::BEST_IN_ROW: {
        best_pos = {stage, best_in_the_row(U, stage)};
        break;
      }
      case SolverMethod::BEST_IN_COLUMN: {
        best_pos = {best_in_the_col(U, stage), stage};
        break;
      }
      case SolverMethod::DO_NOT_TOUCH: {
        best_pos = {stage, stage};
        break;
      }
      default:throw logic_error("TBD");
    }

    swap_rows(U, stage, best_pos.first);
    swap_columns(U, stage, best_pos.second);

//    if (stage != best_pos.first) {
//      cout << "change rows fact " << stage << " " << best_pos.first << "\n";
//    }
//    if (stage != best_pos.second) {
//      cout << "change column fact " << stage << " " << best_pos.second << "\n";
//    }

    Matrix<T> P(I), D(I);
    swap(P[stage], P[best_pos.second]);
    swap(D[stage], D[best_pos.first]);
    Ps.push_back(P);
    Ds.push_back(D);

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
    T k = U[stage][stage];
    for (int col = stage; col < n; ++col) {
      U[stage][col] /= k;
    }

    for (int row = stage + 1; row < n; ++row) {
      sub_row(U, row, stage, U[row][stage]);
    }

//    cout << "After stage " << stage << endl;
//    PrintMatrix(U, "U" + to_string(stage));
//    PrintMatrix(L, "L" + to_string(stage));
//    PrintMatrix(P, "P" + to_string(stage));
//    PrintMatrix(D, "D" + to_string(stage));
  }

  Matrix<T> P = I;
  for (int i = 0; i < Ps.size(); ++i) {
    P = mult(P, Ps[i]);
  }

  Matrix<T> D = I;
  for (int i = Ds.size() - 1; i >= 0; --i) {
    D = mult(D, Ds[i]);
  }

  // D = GetReversed(D);
  // P = GetReversed(P);

  return tie(D, L, U, P);
}

vector<double> SolveLinearSystemUsingDLUP(Matrix<double> A, vector<double> B) {
  int n = A.size();

  Matrix<double> D, L, U, P;
  tie(D, L, U, P) = DLUP_Decomposition(A, SolverMethod::BEST_IN_COLUMN);
//  PrintMatrix(D, "D");
//  PrintMatrix(L, "L");
//  PrintMatrix(U, "U");
//  PrintMatrix(P, "P");
//  if (!compareMatrixOfDouble(A, mult(mult(GetReversed(D), L), mult(U, GetReversed(P))))) {
//    cerr << "Smth goes wrong" << endl;
//    assert(false);
//  }
  vector<double> B_original(B);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (D[i][j]) {
        B[i] = B_original[j];
      }
    }
  }
  // L Y = B
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      B[i] -= L[i][j] * B[j];
    }
    B[i] /= L[i][i];
  }

  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j) {
      B[i] -= U[i][j] * B[j];
    }
  }

  return B;
}
 */