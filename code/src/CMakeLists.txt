cmake_minimum_required( VERSION 3.5 )
project( UQCreator )
include( ExternalProject )
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -frename-registers -DNDEBUG" )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall" )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin/" )

set( UQCreator_SRCS
    main.cpp
    problem.cpp
    burgers.cpp
    polynomial.cpp
    legendre.cpp
    hermite.cpp
    mathtools.cpp
    mesh.cpp
    closure.cpp
    momentsolver.cpp
    limiter.cpp
    nolimiter.cpp
    minmod.cpp
    timesolver.cpp
    thetamethod.cpp
    plotengine.cpp
    gnuplotlib.cpp
    matplotlib.cpp
    sspmultistep.cpp
    heun.cpp
    noplot.cpp
)

add_executable( UQCreator ${UQCreator_SRCS} )

find_package( OpenMP REQUIRED )
if( OPENMP_FOUND )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

find_package( Boost COMPONENTS system iostreams filesystem REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

find_package( BLAS REQUIRED )
include_directories( ${BLAS_INCLUDE_DIR} )

find_package( LAPACK REQUIRED )
include_directories( ${LAPACK_INCLUDE_DIR} )

find_package( PythonLibs 2.7 )
target_include_directories( UQCreator PRIVATE ${PYTHON_INCLUDE_DIRS} )

find_package( Git REQUIRED )
ExternalProject_Add( blaze
  PREFIX ${CMAKE_BINARY_DIR}/blaze
  GIT_REPOSITORY https://bitbucket.org/blaze-lib/blaze
  UPDATE_COMMAND git pull
  UPDATE_DISCONNECTED 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND cmake ../blaze -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/blaze
  INSTALL_COMMAND make install
)
include_directories( ${CMAKE_BINARY_DIR}/blaze/include )
add_dependencies ( UQCreator blaze )

ExternalProject_Add( cpptoml
  PREFIX ${CMAKE_BINARY_DIR}/cpptoml
  GIT_REPOSITORY https://github.com/skystrife/cpptoml
  UPDATE_COMMAND git pull
  UPDATE_DISCONNECTED 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)
include_directories( ${CMAKE_BINARY_DIR}/cpptoml/src/cpptoml/include )
add_dependencies ( UQCreator cpptoml )

ExternalProject_Add( gnuplot-iostream
  PREFIX ${CMAKE_BINARY_DIR}/gnuplot-iostream
  GIT_REPOSITORY https://github.com/dstahlke/gnuplot-iostream
  UPDATE_COMMAND git pull
  UPDATE_DISCONNECTED 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)
include_directories( ${CMAKE_BINARY_DIR}/gnuplot-iostream/src/gnuplot-iostream )
add_dependencies ( UQCreator gnuplot-iostream )

ExternalProject_Add( matplotlib-cpp
  PREFIX ${CMAKE_BINARY_DIR}/matplotlib-cpp
  GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
  UPDATE_COMMAND git pull
  UPDATE_DISCONNECTED 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)
include_directories( ${CMAKE_BINARY_DIR}/matplotlib-cpp/src/matplotlib-cpp )
add_dependencies ( UQCreator matplotlib-cpp )

target_link_libraries( UQCreator LINK_PUBLIC ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${PYTHON_LIBRARIES} )
