include_directories (.)

find_package (Boost COMPONENTS
	graph
	filesystem
	unit_test_framework
	REQUIRED)

link_directories (${Boost_LIBRARY_DIRS})
include_directories (${Boost_INCLUDE_DIRS})

file (GLOB QASSEMBLER_TESTS_SRC *Test.cc)

foreach (test ${QASSEMBLER_TESTS_SRC})
	get_filename_component (basename ${test} NAME_WE)
	add_executable (${basename}.bin ${test})
	target_link_libraries (${basename}.bin ${Boost_LIBRARIES} qassembler)
endforeach(test)

# add_executable (test-qassembler.bin QAssemblerTestSuite.cc ${QASSEMBLER_TESTS_SRC})
# target_link_libraries (test-qassembler.bin ${Boost_LIBRARIES} qassembler)

enable_testing ()
foreach (test ${QASSEMBLER_TESTS_SRC})
	get_filename_component (basename ${test} NAME_WE)
	add_test (${basename} ${basename}.bin)
endforeach(test)
