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

include(ExternalProject)
ExternalProject_Add(blaze
  GIT_REPOSITORY https://bitbucket.org/blaze-lib/blaze
  UPDATE_COMMAND git pull
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/blaze
)

include_directories(${CMAKE_BINARY_DIR}/blaze/include)
link_directories(${CMAKE_BINARY_DIR}/blaze/share)
add_executable ( UQCreator ${UQCreator_SRCS} )
add_dependencies ( UQCreator blaze )
