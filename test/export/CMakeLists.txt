cmake_minimum_required(VERSION 3.18)
project(test-test Fortran)

find_package(dftd4 CONFIG REQUIRED)

add_executable(test-test "example.f90")
target_link_libraries(test-test PUBLIC dftd4::dftd4)
