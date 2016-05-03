cmake_minimum_required (VERSION 2.6)
INCLUDE (CheckIncludeFiles)
include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/AppendCompilerFlags.cmake")

project (Parallel-Range-Light CXX)
set(PROJECT_VENDOR "Julian Labeit")
set(PROJECT_CONTACT "julianlabeit@gmail.com")
set(PROJECT_URL "https://github.com/jlabeit/parallel-range-light")
set(PROJECT_DESCRIPTION "Parallel-range-light: Parallel suffix array algorithm.")
set(CMAKE_BUILD_TYPE "Release")

# Settings about parallization.
option (CILK "Use intel implementation of cilkplus for parallization." OFF) 
option (CILKP "Use gcc implementation of cilkplus for parallization." OFF) 
option (OPENMP "Use OpenMP for parallization for parallization." OFF) 

# Check settings.
if ((CILKP AND OPENMP) OR (CILKP AND CILK) OR (OPENMP AND CILK))
	message(WARNING "Multiple parallelization techniques are set. This may result in problems.")
endif ((CILKP AND OPENMP) OR (CILKP AND CILK) OR (OPENMP AND CILK))
if (NOT(CILK OR CILKP OR OPENMP))
	message(WARNING "No parallelization technique is set, code will be executed sequentially.")
endif (NOT(CILK OR CILKP OR OPENMP))

# Settings for GCC.
if (CILKP)
	append_cxx_compiler_flags("-fcilkplus" "GCC" CMAKE_CXX_FLAGS)
endif (CILKP)
if (OPENMP)
	append_cxx_compiler_flags("-fopenmp" "GCC" CMAKE_CXX_FLAGS)
endif (OPENMP)
if (CILK)
	message( SEND_ERROR "Intel CilkPlus setting is not supported by GCC.")
endif (CILK)
append_cxx_compiler_flags("-Wall -std=c++11" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O2" "GCC" CMAKE_CXX_OPT_FLAGS)

add_subdirectory (include)
add_subdirectory (lib)
add_subdirectory (demo)
add_subdirectory (tests)
