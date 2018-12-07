cmake_minimum_required( VERSION 3.5 )
project( UQCreator )
set( CMAKE_CXX_STANDARD 17 )

set( TAUPATH /home/jannick/Downloads/tau-2.28/x86_64/bin )

set ( CMAKE_CXX_COMPILER ${TAUPATH}/tau_cxx.sh )
set ( CMAKE_C_COMPILER ${TAUPATH}/tau_cc.sh )

set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG" )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -march=native -pg -no-pie" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin/" )

set( SRCS
 ./quadrangle.cpp
 ./matrix.cpp
 ./euler.cpp
 ./burgers.cpp
 ./line.cpp
 ./problem.cpp
 ./legendre.cpp
 ./eulerclosure2d.cpp
 ./cell.cpp
 ./triangle.cpp
 ./mathtools.cpp
 ./euler2d.cpp
 ./timesolver.cpp
 ./closure.cpp
 ./vector.cpp
 ./main.cpp
 ./boundedbarrier.cpp
 ./momentsolver.cpp
 ./mesh.cpp
 ./mesh1d.cpp
 ./stochasticgalerkin.cpp
 ./mesh2d.cpp
 ./expliciteuler.cpp
 ./settings.cpp
 ./polynomial.cpp
 ./vectorspace.cpp
 ./hermite.cpp
 ./eulerclosure.cpp
)

add_executable( ${CMAKE_PROJECT_NAME} ${SRCS} )

find_package( OpenMP REQUIRED )
if( OPENMP_FOUND )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

find_package( BLAS REQUIRED )
include_directories( ${BLAS_INCLUDE_DIR} )

find_package( LAPACK REQUIRED )
include_directories( ${LAPACK_INCLUDE_DIR} )

find_package( VTK REQUIRED COMPONENTS vtkIOXML vtkFiltersCore )

include_directories( ${CMAKE_SOURCE_DIR}/../ext/cpptoml/include )
include_directories( ${CMAKE_SOURCE_DIR}/../ext/spdlog/include )

target_link_libraries( ${CMAKE_PROJECT_NAME} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${VTK_LIBRARIES} -lstdc++fs )
