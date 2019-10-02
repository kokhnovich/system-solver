//
// Created by user on 25.09.19.
//

/*
 * task1
 * task2
 * task3... Done!
 * task4
 * task5
 */

void SolveTask1() {
  auto solver = new HW_Solver<double>();
  Matrix<double> A1 = {
      {-3, 1, 0, 0},
      {0, 2, 5, 0},
      {-3, -4, -5, -2},
      {3, -2, 1, -1}
  };

  Matrix<double> A2 = {
      {-5, -5, 0, 0, 0, 0, 0, 0, 0},
      {3, 3, 0, 0, 0, 0, 0, 0, 0},
      {-5, 0, -4, 5, 0, 0, 0, 0, 0},
      {-4, 2, -5, 4, 0, 0, 0, 0, 0},
      {-3, -4, 1, 4, -4, 3, 0, 0, 0},
      {-5, 3, 4, 2, 4, -2, 2, 0, 0},
      {4, 0, -1, -1, 3, 1, 6, -5, 0},
      {0, 2, 2, -5, -2, 0, 1, 3, -4},
      {-5, 3, -5, 3, 1, 0, -3, 1, -4}
  };

//  try {
//    auto revA2 = solver->task2_gauss(A2);
//    PrintMatrix(revA2);
//  } catch (exception& e) {
//    cerr << e.what() << endl;
//  }

//  auto revA1 = solver->task2_gauss(A1);
//  if (!compareMatrixOfDouble(mult(A1, revA1), getIdentityMatrix<double>(A1.size()))) {
//    PrintMatrix(A1, "A1 original");
//    PrintMatrix(mult(A1, revA1), "A * A^{-1} gmo");
//  }

  vector<int> times;
  for (int cnt = 1800; cnt <= 2300; cnt += 200) {
    Matrix<double> A(solver->task2_random_strange_matrix(cnt));
    Matrix<double> B(getIdentityMatrix<double>(cnt));
    auto timer = new LogDuration("For size=" + to_string(cnt) + " is ");
    Matrix<double> ans = solver->task2_gauss(A);
    times.push_back(timer->getTime());
    delete timer;
//    if (!compareMatrixOfDouble(mult(A, ans), B)) {
//      // PrintPartOfMatrix(mult(A, ans), 10, "A * A^{-1}");
//      // PrintPartOfMatrix(ans, 10, "ans");
//      cerr << "wroong" << endl;
//    }
  }

  for (const auto& time : times) {
    cout << time << " ";
  }
}

void SolveTask2() {
  cout << "Solving task 2\n";
  Matrix<double> A = {
      {-3, 4, 2, -3, 2, 3, 4, 3},
      {-3, 4, -1, 1, -4, -5, -1, 0},
      {-6, 8, 1, -1, 1, 3, -3, -5},
      {-12, 16, 2, -3, 4, -2, -3, 4},
      {-24, 32, 4, -6, 3, -3, -1, -2},
      {-48, 64, 8, -12, 6, -4, 3, -3},
      {-96, 128, 16, -24, 12, -8, -1, 4},
      {2, -1, -3, -4, -5, 4, 1, -5}
  };
  vector<double> b = {79, -51, -29, 33, 2, 59, 149, -59};
  // Matrix<double> B = make2Dfrom1D(b);
  auto solver = new Solver<double>();
  PrintMatrix(solver->SolveLinearSystemUsingDLUP(A, b), "Ax=b");
  // 1 2 3 4 5 6 7 8

  Matrix<double> A2 = {
      {1, 1, 1, 1, 1, 1, 1, 1, 1},
      {1, 2, 4, 8, 16, 32, 64, 128, 256},
      {1, 3, 9, 27, 81, 243, 729, 2187, 6561},
      {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536},
      {1, 5, 25, 125, 625, 3125, 15625, 78125, 390625},
      {1, 6, 36, 216, 1296, 7776, 46656, 279936, 1679616},
      {1, 7, 49, 343, 2401, 16807, 117649, 823543, 5764801},
      {1, 8, 64, 512, 4096, 32768, 262144, 2097152, 16777216},
      {1, 9, 81, 729, 6561, 59049, 531441, 4782969, 43046721}
  };
  vector<double> b2 = {9, 511, 9841, 87381, 488281, 2015539, 6725601, 19173961, 48427561};

  PrintMatrix(solver->SolveLinearSystemUsingDLUP(A2, b2), "A2 x = b2");
}

void SolveTask3() {
  cout << "Solving task 3\n";
  auto solver = new Solver<double, greater_using_abs<double>>();
  std::uniform_int_distribution<std::mt19937_64::result_type> udist(-5, 5);
  std::mt19937_64 generator(std::random_device{}());

//  Matrix<double> A(generateRandomSymmetricMatrix<double>(10, 4, 14));
//  vector<double> x(10, 1.);
//  vector<double> B(generateAnsMatrix(A, x));
//  PrintMatrix(solver->SolveLinearSystemUsingLDLt(A, B), "ldlt x = b");
//  PrintMatrix(solver->SolveLinearSystemUsingDLUP(A, B), "dlup x = b");
//  PrintMatrix(A, "A");
//  return;

  vector<int> times;
  for (int cnt = 100; cnt <= 2000; cnt += 100) {
    Matrix<double> A(generateRandomSymmetricMatrix<double>(cnt, 2, 10));
    // PrintMatrix(A, "A");
    vector<double> x(cnt, 1.);
    vector<double> B(generateAnsMatrix(A, x));
    auto timer = new LogDuration("For size=" + to_string(cnt) + " is ");
    auto ans = solver->SolveLinearSystemUsingLDLt(A, B);
    times.push_back(timer->getTime());
    delete timer;
    if (!compareVectorOfDouble(x, ans)) {
      PrintMatrix(x, "x");
      PrintMatrix(ans, "ans");
    }
  }
  for (auto& i : times) {
    cout << i << " ";
  }
}