CMAKE_MINIMUM_REQUIRED(VERSION 3.24)

PROJECT(Bernoulli_Euler_beam)

SET( CMAKE_CXX_COMPILER "g++")

ADD_COMPILE_OPTIONS(-std=c++20)

SET(Eigen3_DIR CMAKE_INSTALL_PREFIX/share/eigen3/cmake)

FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

INCLUDE_DIRECTORIES("/usr/include/eigen3")

ADD_EXECUTABLE(Bernoulli_Euler_beam main.cpp Tools.cpp)

TARGET_LINK_LIBRARIES(Bernoulli_Euler_beam Eigen3::Eigen)

