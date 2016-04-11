# GTest.
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../lib")

add_executable(runTests tests.cpp)
target_link_libraries (runTests libprange divsufsort cilkrts ${GTEST_LIBRARIES} pthread)
