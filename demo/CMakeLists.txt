add_executable(PRangeLightDemo main.cpp)
target_link_libraries (PRangeLightDemo libprange cilkrts)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include")
