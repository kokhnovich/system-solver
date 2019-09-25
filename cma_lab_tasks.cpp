//
// Created by user on 25.09.19.
//

void SolveTask1() {

}

void SolveTask2() {

}

void SolveTask3() {
  auto solver = new Solver<double>();
  std::uniform_int_distribution<std::mt19937_64::result_type> udist(-5, 5);
  std::mt19937_64 generator(std::random_device{}());
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