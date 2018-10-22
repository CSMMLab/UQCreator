cmake_minimum_required( VERSION 3.5 )
project( UQCreator )
include( ExternalProject )
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG" )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -march=native -pg -no-pie" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall" )
add_definitions( -D_GLIBCXX_PARALLEL -DQCUSTOMPLOT_USE_OPENGL )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin/" )

set( CMAKE_AUTOMOC ON )
set( CMAKE_AUTOUIC ON )
set( CMAKE_AUTORCC ON )
set( CMAKE_INCLUDE_CURRENT_DIR ON )

file(GLOB_RECURSE UQCreator_SRCS RELATIVE ${CMAKE_SOURCE_DIR} "*.cpp")
add_executable( ${CMAKE_PROJECT_NAME} ${UQCreator_SRCS} )

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

find_package( VTK REQUIRED COMPONENTS vtkIOXML vtkFiltersCore )

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
add_dependencies ( ${CMAKE_PROJECT_NAME} blaze )

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
add_dependencies ( ${CMAKE_PROJECT_NAME} cpptoml )

find_package( Qt5 COMPONENTS Core Widgets Gui PrintSupport REQUIRED )

set( OpenGL_GL_PREFERENCE LEGACY )
find_package( OpenGL )

target_link_libraries( ${CMAKE_PROJECT_NAME}  LINK_PUBLIC Qt5::Core Qt5::Widgets Qt5::Gui Qt5::PrintSupport ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${VTK_LIBRARIES} -lstdc++fs)
