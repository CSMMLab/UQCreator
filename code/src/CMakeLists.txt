cmake_minimum_required( VERSION 3.7 )
project( UQCreator )
set ( CMAKE_BUILD_TYPE Release )
add_definitions ( -std=c++14 -Wall -O3 )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")

set ( UQCreator_HDR
    )
set ( UQCreator_SRCS
	  main.cpp
	)

add_executable ( UQCreator ${UQCreator_SRCS} )

include(ExternalProject)
ExternalProject_Add(blaze
  GIT_REPOSITORY https://bitbucket.org/blaze-lib/blaze
  UPDATE_COMMAND git pull
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/blaze
)
include_directories(${CMAKE_BINARY_DIR}/blaze/include)
add_dependencies ( UQCreator blaze )

ExternalProject_Add(cpptoml
  GIT_REPOSITORY https://github.com/skystrife/cpptoml.git 
  UPDATE_COMMAND git pull
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/cpptoml
)
include_directories(${CMAKE_BINARY_DIR}/cpptoml/include)
add_dependencies ( UQCreator cpptoml )
