cmake_minimum_required(VERSION 3.17)
project(BcProject C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99 -Wall -pedantic")

add_executable(BcProject main.c helper_functions.c helper_functions.h fitness_functions.c fitness_functions.h parameters.c parameters.h search.c search.h)

target_link_libraries(BcProject PRIVATE m)
