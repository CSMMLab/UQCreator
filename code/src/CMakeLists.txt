cmake_minimum_required( VERSION 3.8 )
project( UQCreator )
set ( CMAKE_BUILD_TYPE Release )
set (CMAKE_CXX_STANDARD 14)
add_definitions ( -std=c++14 -Wall -ansi -O3 -DNDEBUG -mavx )

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
          polynomial.cpp
          closure.cpp
          momentsolver.cpp
          timediscretization.cpp
	)

add_executable ( UQCreator ${UQCreator_SRCS})

find_package( Boost COMPONENTS system iostreams REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include(ExternalProject)
ExternalProject_Add(blaze
  PREFIX ${CMAKE_BINARY_DIR}/blaze
  GIT_REPOSITORY https://bitbucket.org/blaze-lib/blaze
  UPDATE_COMMAND git pull
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)
include_directories(${CMAKE_BINARY_DIR}/blaze/src/blaze)
add_dependencies ( UQCreator blaze )

find_package( LAPACK REQUIRED )
include_directories( ${LAPACK_INCLUDE_DIR} )

ExternalProject_Add(cpptoml
  PREFIX ${CMAKE_BINARY_DIR}/cpptoml
  GIT_REPOSITORY https://github.com/skystrife/cpptoml
  UPDATE_COMMAND git pull
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)
include_directories(${CMAKE_BINARY_DIR}/cpptoml/src/cpptoml/include)
add_dependencies ( UQCreator cpptoml )

TARGET_LINK_LIBRARIES( UQCreator LINK_PUBLIC ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
