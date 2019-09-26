//
// Created by user on 25.09.19.
//

void SolveTask1() {

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
  auto solver = new Solver<double>();
  std::uniform_int_distribution<std::mt19937_64::result_type> udist(-5, 5);
  std::mt19937_64 generator(std::random_device{}());

  Matrix<double> A(generateRandomSymmetricMatrix<double>(10, 4, 14));
  // PrintMatrix(A, "A");
  vector<double> x(10, 1.);
  vector<double> B(generateAnsMatrix(A, x));
  solver->SolveSystemUsingLDLt(A, make2Dfrom1D(B));
  return;

  vector<int> times;
  for (int cnt = 100; cnt <= 2000; cnt += 100) {
    Matrix<double> A(generateRandomSymmetricMatrix<double>(cnt, 0, 10));
    // PrintMatrix(A, "A");
    vector<double> x(cnt, 1.);
    vector<double> B(generateAnsMatrix(A, x));
    auto timer = new LogDuration("For size=" + to_string(cnt) + " is ");
    auto ans = solver->SolveSystemUsingLU(A, make2Dfrom1D(B), SolverMethod::BEST_IN_ROW);
    times.push_back(timer->getTime());
    delete timer;
    if (!compareVectorOfDouble(x, make1Dfrom2D(ans))) {
      PrintMatrix(x, "x");
      PrintMatrix(make1Dfrom2D(ans), "ans");
    }
  }
  for (auto& i : times) {
    cout << i << " ";
  }
}