cmake_minimum_required( VERSION 3.5 )
project( UQCreator )
set( CMAKE_CXX_STANDARD 17 )

set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native" )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -pg -no-pie" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin/" )

file(GLOB_RECURSE SRCS RELATIVE ${CMAKE_SOURCE_DIR} "*.cpp")
add_executable( ${CMAKE_PROJECT_NAME} ${SRCS} )

find_package( OpenMP REQUIRED )
if( OPENMP_FOUND )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )

find_package( LAPACK REQUIRED )
include_directories( ${LAPACK_INCLUDE_DIR} )

find_package( VTK REQUIRED COMPONENTS vtkIOGeometry vtkFiltersCore )

include_directories( ${CMAKE_SOURCE_DIR}/../ext/cpptoml/include )
include_directories( ${CMAKE_SOURCE_DIR}/../ext/spdlog/include )

execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions('-DGIT_HASH="${GIT_HASH}"')

target_link_libraries( ${CMAKE_PROJECT_NAME} ${LAPACK_LIBRARIES} ${VTK_LIBRARIES} ${MPI_LIBRARIES} -lstdc++fs)
