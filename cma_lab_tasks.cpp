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
    auto timer = new LogDuration("For size=" + to_string(cnt) + " is ");
    Matrix<double> A(generateRandomSymmetricMatrix<double>(cnt, 0, 10));
    // PrintMatrix(A, "A");
    vector<double> x(cnt, udist(generator));
    vector<double> B(generateAnsMatrix(A, x));
    auto ans = solver->SolveSystemUsingLU(A, make2Dfrom1D(B));
    // PrintMatrix(ans, "for " + to_string(cnt));
    times.push_back(timer->getTime());
    delete timer;
  }
  for(auto& i : times) {
    cout << i << " ";
  }
}