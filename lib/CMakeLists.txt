add_library(libprange parallel-range.cpp)
install(TARGETS libprange
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)
