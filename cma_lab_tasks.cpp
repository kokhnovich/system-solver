//
// Created by user on 25.09.19.
//

void SolveTask1() {

}

void SolveTask2() {

}


void SolveTask3() {
  auto solver = new Solver<Fraction>();
  for(int cnt = 100; cnt <= 2000; ++cnt) {
    auto timer = new LogDuration("For size="+to_string(cnt)+" is ");
    // solver->SolveSystem();
    delete timer;
  }
}