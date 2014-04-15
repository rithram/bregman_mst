cmake_minimum_required(VERSION 2.6)

project(BregmanMST)

find_package(Boost 
  COMPONENTS
    program_options
  REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x -O0 -g -ggdb")
#if(DEBUG)
#  add_definitions(-DDEBUG)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -ggdb -Wall")
#else()
#  add_definitions(-DNDEBUG)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wextra")
#endif(DEBUG)

add_executable(test_data_main 
  test_data_main.cpp)

add_executable(test_divergences 
  test_divergences.cpp)

add_executable(test_bregman_ball 
  test_bregman_ball.cpp)

add_executable(test_kmeans_splitter 
  test_kmeans_splitter.cpp)

add_executable(test_bbtree_main 
  test_bbtree_main.cpp)

add_executable(test_leftnn_search 
  test_leftnn_search.cpp)

add_executable(test_search_main 
  test_search_main.cpp)
target_link_libraries(test_search_main 
  ${Boost_LIBRARIES})

#add_executable(test_mst 
#  test_mst.cpp  union_find.cpp)
#target_link_libraries(test_mst
#  ${Boost_LIBRARIES})
