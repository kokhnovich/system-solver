cmake_minimum_required(VERSION 3.14)
project(system-solver)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS -pthread)

add_executable(system-solver main.cpp)
