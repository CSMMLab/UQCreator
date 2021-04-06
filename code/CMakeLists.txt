cmake_minimum_required( VERSION 3.16 )
project( UQCreator VERSION 1.0.0 LANGUAGES CXX )


### OPTIONS #####################################
option( UQCREATOR_BUILD_SHARED_LIB "Enables building of shared library for usage with e.g. the Python interface." OFF )
option( UQCREATOR_BUILD_UNITY "Enables unity build for faster compile times." OFF )
#################################################


### COMPILER ####################################
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( UQCREATOR_RELEASE_OPTIONS -march=native -w )
set( UQCREATOR_RELWITHDEBINFO_OPTIONS -march=native -pg -no-pie )
set( UQCREATOR_DEBUG_OPTIONS -Wall -Wextra -Wpedantic )
if( BUILD_UNITY )
    message( STATUS "Unity build enabled" )
    set( CMAKE_UNITY_BUILD ON )
    set( CMAKE_UNITY_BUILD_BATCH_SIZE 0 )
endif()
#################################################


### LIBRARIES ###################################
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
#################################################


### EXTERNAL LIBRARIES ##########################
include_directories( ${CMAKE_SOURCE_DIR}/ext/cpptoml/include )
include_directories( ${CMAKE_SOURCE_DIR}/ext/spdlog/include )
#################################################


### MISC ########################################
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_compile_definitions( GIT_HASH="${GIT_HASH}" )
#################################################


### BUILD UQCREATOR #############################
file( GLOB_RECURSE SRCS RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp" "include/*.h" )
include_directories( ${CMAKE_SOURCE_DIR}/include )
if( UQCREATOR_BUILD_SHARED_LIB )
    add_library( ${CMAKE_PROJECT_NAME} SHARED ${SRCS} )
else()
    add_executable( ${CMAKE_PROJECT_NAME} ${SRCS} )
endif()
target_link_libraries( ${CMAKE_PROJECT_NAME} ${LAPACK_LIBRARIES} ${VTK_LIBRARIES} ${MPI_LIBRARIES} -lstdc++fs )
set_target_properties( ${CMAKE_PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin" )
target_compile_options( ${CMAKE_PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:${UQCREATOR_DEBUG_OPTIONS}>" )
target_compile_options( ${CMAKE_PROJECT_NAME} PUBLIC "$<$<CONFIG:RELWITHDEBINFO>:${UQCREATOR_RELWITHDEBINFO_OPTIONS}>" )
target_compile_options( ${CMAKE_PROJECT_NAME} PUBLIC "$<$<CONFIG:RELEASE>:${UQCREATOR_RELEASE_OPTIONS}>" )
#################################################
