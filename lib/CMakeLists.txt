include_directories (.)

find_package (Boost COMPONENTS 
	filesystem
	graph
	REQUIRED)
link_directories (${Boost_LIBRARY_DIRS})
include_directories (${Boost_INCLUDE_DIRS})

find_package (ZLIB REQUIRED)
link_directories (${ZLIB_LIBRARY_DIRS})
include_directories (${ZLIB_INCLUDE_DIRS})

file (GLOB_RECURSE QASSEMBLER_SRC *.cc *.hh)

add_library (qassembler SHARED ${QASSEMBLER_SRC})
target_link_libraries (qassembler ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${EXTRA_LINK_LIBRARIES})

install (TARGETS qassembler LIBRARY DESTINATION lib)
