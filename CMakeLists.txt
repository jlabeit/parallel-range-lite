cmake_minimum_required (VERSION 2.6)
project (Parallel-Range-Light)
add_executable(Parallel-Range-Light main.cpp)

#Compiler Flags.
target_link_libraries(Parallel-Range-Light divsufsort cilkrts)
add_definitions(-g -Wall -std=c++11 -fcilkplus)
#add_definitions(-Wall -std=c++11 -DCILKP -fcilkplus -lcilkplus -lcilkrts)

# GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(runTests tests.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)
