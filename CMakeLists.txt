cmake_minimum_required (VERSION 2.6)
include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/AppendCompilerFlags.cmake")

project (Parallel-Range-Light CXX)
set(PROJECT_VENDOR "Julian Labeit")
set(PROJECT_CONTACT "julianlabeit@gmail.com")
set(PROJECT_URL "https://github.com/jlabeit/parallel-range-light")
set(PROJECT_DESCRIPTION "Parallel-range-light: Parallel suffix array algorithm.")
set(CMAKE_BUILD_TYPE "Release")

# Settings about parallization.
option (CILK "Use intel implementation of cilkplus for parallization." OFF) 
option (CILKP "Use gcc implementation of cilkplus for parallization." ON) 
option (OPENMP "Use OpenMP for parallization for parallization." OFF) 

append_cxx_compiler_flags("-Wall -std=c++11 -fcilkplus" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O2" "GCC" CMAKE_CXX_OPT_FLAGS)

add_subdirectory (include)
add_subdirectory (lib)
add_subdirectory (demo)
add_subdirectory (tests)
