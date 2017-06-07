cmake_minimum_required( VERSION 3.5 )
project( UQCreator )
set ( CMAKE_BUILD_TYPE Release )
add_definitions ( -std=c++14 -Wall -O3 )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")

set ( UQCreator_HDR
    )
set ( UQCreator_SRCS
	  main.cpp
	)

find_package( blaze QUIET )
if( blaze_FOUND )
    add_library( blaze_target INTERFACE )
    target_link_libraries( blaze_target INTERFACE blaze::blaze )
    add_executable ( UQCreator ${UQCreator_SRCS} )
    target_link_libraries( UQCreator blaze_target )
else()
    include( ${CMAKE_CURRENT_SOURCE_DIR}/../lib/blaze-3.1/CMakeLists.txt )
endif()

