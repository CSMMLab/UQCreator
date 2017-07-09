cmake_minimum_required( VERSION 3.7 )
project( UQCreator )
set ( CMAKE_BUILD_TYPE Release )
set (CMAKE_CXX_STANDARD 14)
add_definitions ( -std=c++14 -Wall -O3 )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")

set ( UQCreator_SRCS
	  main.cpp
          problem.cpp
          quadrature.cpp
          legendre.cpp
          hermite.cpp
          mathtools.cpp
          mesh.cpp
          burgers.cpp
          euler.cpp
	)

add_executable ( UQCreator ${UQCreator_SRCS})

find_package( Boost COMPONENTS system iostreams REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

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

TARGET_LINK_LIBRARIES( UQCreator LINK_PUBLIC ${Boost_LIBRARIES} )
