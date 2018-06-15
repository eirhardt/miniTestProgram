cmake_minimum_required(VERSION 3.10)
project(miniTestProgram)

set(CMAKE_CXX_STANDARD 17)

include_directories(basic fem posit utils src)  # try this--works


add_executable(miniTestProgram main.cpp)