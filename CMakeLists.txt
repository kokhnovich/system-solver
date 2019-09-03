cmake_minimum_required(VERSION 3.14)
project(system-solver)

set(CMAKE_CXX_STANDARD 17)

add_executable(system-solver main.cpp test_runner.cpp test_runner.h)
