include_directories (.)

find_package (Boost COMPONENTS 
	system
	filesystem
	program_options
	REQUIRED)
link_directories (${Boost_LIBRARY_DIRS})
include_directories (${Boost_INCLUDE_DIRS})

add_executable (denovo-qassembler QAssembler.cc)
target_link_libraries (denovo-qassembler ${Boost_LIBRARIES} qassembler ${EXTRA_LINK_LIBRARIES})

install (PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/denovo-qassembler" DESTINATION bin)
