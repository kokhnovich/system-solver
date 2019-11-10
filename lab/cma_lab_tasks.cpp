//
// Created by user on 25.09.19.
//

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

  auto revA1 = solver->task1_gauss(A1);
  if (!compareMatrixOfDouble(mult(A1, revA1), getIdentityMatrix<double>(A1.size()))) {
    PrintMatrix(A1, "A1 original");
    PrintMatrix(mult(A1, revA1), "A * A^{-1} my");
  } else {
    cout << "1. OK" << endl;
  }

  /* should't work because of non-existance of solution
  auto revA2 = solver->task1_gauss(A2);
  if (!compareMatrixOfDouble(mult(A2, revA2), getIdentityMatrix<double>(A2.size()))) {
    PrintMatrix(A1, "A2 original");
    PrintMatrix(mult(A1, revA1), "A * A^{-1} my");
  } else {
    cout << "2. OK" << endl;
  }
  */


  vector<int> times;
  for (int cnt = 7000; cnt <= 7500; cnt += 100) {
    Matrix<double> A(solver->task1_random_strange_matrix(cnt));
    Matrix<double> B(getIdentityMatrix<double>(cnt));
    auto timer = new LogDuration("For size=" + to_string(cnt) + " is ");
    try {
      Matrix<double> ans = solver->task1_gauss(A);
      times.push_back(timer->getTime());
      delete timer;
      auto test_timer = new LogDuration("Testing takes ");
      /*
      if (!compareMatrixOfDouble(mult(A, ans), B)) {
        cerr << "ooops... smth goes wrong" << endl;
        PrintMatrix(ans, "ans");
        PrintMatrix(A, "A");
        PrintMatrix(B, "B");
      }
      delete test_timer;
      */
    } catch (std::exception& e) {
      cerr << e.what() << endl;
      continue;
    }
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
  auto solver = new HW_Solver<double>();
  PrintMatrix(A, "A1");
  PrintMatrix(b, "b1");
  PrintMatrix(solver->task2_dlup(A, b), "Ax=b");
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

  PrintMatrix(A2, "A2");
  PrintMatrix(b2, "b2");
  PrintMatrix(solver->task2_dlup(A2, b2), "A2 x = b2");
  // 1 1 1 1 1 1 1 1

  std::uniform_int_distribution<std::mt19937_64::result_type> udist(0, 8);
  std::mt19937_64 generator(std::random_device{}());

  int tests = 25;
  while (tests--) {
    int rnd = udist(generator);
    b2[rnd] *= (rnd % 2 ? 1.0001 : 0.9999);
    PrintMatrix(b2, "b2");
    PrintMatrix(solver->task2_dlup(A2, b2), "test" + to_string(tests));
  }
}

void SolveTask3() {
  cout << "Solving task 3\n";
  auto solver = new HW_Solver<double>();

  /*
  std::uniform_int_distribution<std::mt19937_64::result_type> udist(-5, 5);
  std::mt19937_64 generator(std::random_device{}());
  //Matrix<double> A(generateRandomSymmetricMatrix<double>(4, 2, 10));
  Matrix<double> A = {{1, 2, 3, 4},
                      {2, 5, 8, 9},
                      {3, 8, 6, 10},
                      {4, 9, 10, 7}};
  vector<double> x(4, 1.);
  vector<double> b(generateAnsMatrix(A, x));
  PrintMatrix(b, "b");
  PrintMatrix(A, "A");
  PrintMatrix(solver->task3_ldlt(A, b), "ldlt x = b");
  PrintMatrix(solver->task3_lu_for_sym(A, b), "lu x = b");
  PrintMatrix(A, "A");
  */

  vector<pair<int, int>> times;
  cout << "size --- lu --- ldlt" << endl;
  for (int cnt = 100; cnt <= 2000; cnt += 100) {
    Matrix<double> A(generateRandomSymmetricMatrix<double>(cnt, 2, 100));
    int time1 = 0, time2 = 0;
    vector<double> x(cnt, 3.);
    vector<double> B(generateAnsMatrix(A, x));
    {
      Matrix<double> AA(A);
      vector<double> BB(B);
      auto timer = new LogDuration("LDLt. Size=" + to_string(cnt) + " is ");
      auto ans = solver->task3_ldlt(AA, BB);
      time1 = timer->getTime();
      delete timer;
      if (!compareVectorOfDouble(x, ans)) {
        cerr << "WRONG LDLt" << endl;
        PrintMatrix(x, "x");
        PrintMatrix(ans, "ans");
      }
    }
    {
      auto timer = new LogDuration("LU. Size=" + to_string(cnt) + " is ");
      auto ans = solver->task3_lu_for_sym(A, B);
      time2 = timer->getTime();
      delete timer;
      if (!compareVectorOfDouble(x, ans)) {
        cerr << "WRONG LU" << endl;
        PrintMatrix(x, "x");
        PrintMatrix(ans, "ans");
      }
    }
    cerr << time1 << "   " << time2 << endl;
    times.push_back(make_pair(time1, time2));
  }

  cout << "\n\n";
  for (auto& i : times) {
    cout << i.first << " ";
  }
  cout << endl;
  for (auto& i : times) {
    cout << i.second << " ";
  }
  cout << endl;

}

void SolveTask4() {
  cout << "Solving task 4\n";
  auto solver = new Solver<double>();

  Matrix<double> A_1 = {
      {0, -2, 2},
      {2, 3, 5},
      {-2, 4, 1},
      {-2, 3, 0}
  };
  vector<double> b_1 = {0, 10, 3, 1};
  PrintMatrix(solver->SolveThreeDiagonalSystem(A_1, b_1), "Solution");

  Matrix<double> A_2 = {
      {0, 4, -3},
      {4, -6, 4},
      {2, 4, -4},
      {-2, 3, -1},
      {5, -2, 4},
      {-5, 1, -4},
      {-3, -6, 0},
      {0, 2, 4},
      {0, 3, 2},
      {5, -4, 0}
  };
  vector<double> b_2 = {1, 2, 2, 0, 7, -8, -9, 6, 5, 1};
  PrintMatrix(solver->SolveThreeDiagonalSystem(A_2, b_2), "Solution");
}

void SolveTask5() {
  cout << "Solving task 5\n";
  auto solver = HW_Solver<double>();
  vector<int> ranges = {500, 1000, 2000, 4000};
  for (auto& range : ranges) {
    for (double w = 0.; w <= 2.; w += 0.1) {
      cout << range << " " << fixed << setprecision(2) << w << " ";
      auto ans = solver.task5_relax(range, w);
//      for (int i = 0; i < 10; ++i) {
//        cout << ans.first[i] << " & ";
//      }
//      cout << endl;
      cout << ans.second << endl;
    }
  }
}