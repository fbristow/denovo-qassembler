cmake_minimum_required (VERSION 2.6)
project (QAssembler CXX)

# we use a custom CMake module to find Log4cxx, so add that module to the
# list of cmake modules.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/lib64)

# allow us to use cpack to package the project.
include (CPack)

# version information
set (QAssembler_VERSION_MAJOR 1)
set (QAssembler_VERSION_MINOR 0)

# All packages use Log4cxx, find it here for everyone
find_package (Log4cxx)
if (LOG4CXX_FOUND)
	include_directories (${LOG4CXX_INCLUDE_DIRS})
	link_directories (${LOG4CXX_LIBRARY_DIRS})
	option (USE_LOG4CXX "Use log4cxx for logging?" ON)
	set (EXTRA_LINK_LIBRARIES ${LOG4CXX_LIBRARIES})
endif ()

# optional severity level for logging if log4cxx is not found.
set (LOGGING_SEVERITY "0" CACHE STRING "")

# Both src and tests use log4cxx if it's installed. So we're going to add
# some global configuration stuff into a common config header. 
configure_file (
	"${PROJECT_SOURCE_DIR}/QAssemblerConfig.h.in"
	"${PROJECT_BINARY_DIR}/QAssemblerConfig.h"
)

# Make sure we can find QAssemblerConfig.h when we're building:
include_directories ("${PROJECT_BINARY_DIR}")
# We build most of the objects as a shared library, allow those to be
# found when building.
include_directories (lib)

# Some custom CXXFLAGS options
#add_definitions (-O3 -march=native -mtune=generic -pipe)
add_definitions (-O3 -pipe -Wall)

# add all of the directories for building
add_subdirectory (src)
add_subdirectory (lib)

enable_testing ()
add_subdirectory (tests)

# install (FILES "${PROJECT_BINARY_DIR}/QAssemblerConfig.h" DESTINATION include)

