cmake_minimum_required( VERSION 3.5 )
project( FVSolver1D )
include( ExternalProject )
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "../../lib/" )

set( SRCS
    Vector.cpp
    Cell.cpp
    FVSolver1D.cpp
    Cell1D.cpp
    PhysicalProblem.cpp
    Mesh.cpp
    Mesh1D.cpp
    RiemannSolver.cpp
    HLLSolver1D.cpp
    TimeIntegrator.cpp
    Heun.cpp
    Rhs.cpp
    Rhs1D.cpp
    EulerEquationsIdGas.cpp
    Advection.cpp
    BurgersEquation.cpp
    Limiter.cpp
    VanLeer1D.cpp
    MinMod1D.cpp
    Settings.cpp
)

find_package( Matlab COMPONENTS MX_LIBRARY MEX_COMPILER REQUIRED )
include_directories( ${Matlab_ROOT_DIR}/extern/include/ )

add_library( ${PROJECT_NAME} STATIC ${SRCS} )
add_definitions( -DMATLAB_MEX_FILE )
set_target_properties( ${PROJECT_NAME} PROPERTIES SUFFIX .${Matlab_MEX_EXTENSION})
target_link_libraries( ${PROJECT_NAME} ${MATLAB_LIBRARIES})

