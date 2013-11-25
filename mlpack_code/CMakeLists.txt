cmake_minimum_required(VERSION 2.6)

project(BregmanMST)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -ggdb -std=gnu++0x")




add_executable(test_data_main test_data_main.cpp)


add_executable(test_divergences test_divergences.cpp)


add_executable(test_bregman_ball test_bregman_ball.cpp)

