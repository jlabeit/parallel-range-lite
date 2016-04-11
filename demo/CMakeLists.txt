add_executable(PRangeLightDemo main.cpp)
target_link_libraries (PRangeLightDemo libprange divsufsort cilkrts)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include")
