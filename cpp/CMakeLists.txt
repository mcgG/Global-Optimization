project( cgal )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

find_package(CGAL QUIET COMPONENTS Core )
find_package(GLUT)
find_package(OpenGL)
SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL -framework GLUT")

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include( CGAL_VersionUtils )
  find_package( MPFI )
   if( MPFI_FOUND )
    include( ${MPFI_USE_FILE} )
  endif( MPFI_FOUND )

  include_directories (BEFORE)

  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    include_directories (BEFORE "../../include")

    create_single_source_cgal_program("Delaunay.cpp")

  else()
    message(STATUS "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled.")
  endif()

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
