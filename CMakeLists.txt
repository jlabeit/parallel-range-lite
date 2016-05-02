cmake_minimum_required (VERSION 2.6)
include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/AppendCompilerFlags.cmake")

project (Parallel-Range-Light CXX)
set(PROJECT_VENDOR "Julian Labeit")
set(PROJECT_CONTACT "julianlabeit@gmail.com")
set(PROJECT_URL "https://github.com/jlabeit/parallel-range-light")
set(PROJECT_DESCRIPTION "Parallel-range-light: Parallel suffix array algorithm.")
set(CMAKE_BUILD_TYPE "Release")

add_definitions(-g -O2 -Wall -std=c++11 -fcilkplus -DCILKP)
append_cxx_compiler_flags("-Wall -std=c++11 -fcilkplus -DCILKP" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O2" "GCC" CMAKE_CXX_OPT_FLAGS)

add_subdirectory (include)
add_subdirectory (lib)
add_subdirectory (demo)
add_subdirectory (tests)
