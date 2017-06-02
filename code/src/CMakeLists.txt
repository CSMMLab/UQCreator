project( UQCreator )
cmake_minimum_required( VERSION 2.8 )
set ( CMAKE_BUILD_TYPE Release )
add_definitions ( -std=c++11 -Wall -O3 )
include_directories (
	${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	.)
set ( UQCreator_HDRS
    )
set ( UQCreator_SRCS
	main.cpp
	)

add_executable ( ../bin/UQCreator ${UQCreator_SRCS} )
target_link_libraries ( UQCreator )
