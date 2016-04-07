cmake_minimum_required (VERSION 2.6)
project (Parallel-Ranke-Light)
add_executable(Parallel-Ranke-Light main.cpp)

#Compiler Flags.
add_definitions(-Wall -std=c++11 -fcilkplus -lcilkplus -lcilkrts)
#add_definitions(-Wall -std=c++11 -DCILKP -fcilkplus -lcilkplus -lcilkrts)

# GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(runTests tests.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)
