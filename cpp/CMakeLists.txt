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

    create_single_source_cgal_program( "demo.cpp" )
    create_single_source_cgal_program( "demo1.cpp" )
    create_single_source_cgal_program( "voronoi.cpp" )
    create_single_source_cgal_program( "conforming.cpp" )
    create_single_source_cgal_program( "mesh.cpp" )
    create_single_source_cgal_program( "demo2.cpp" )
    create_single_source_cgal_program("Global_Optimization.cpp")
    create_single_source_cgal_program("Delaunay.cpp")
    create_single_source_cgal_program("D3.cpp")
    create_single_source_cgal_program("Delaunay_d.cpp")
  else()
    message(STATUS "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled.")
  endif()

  target_link_libraries(demo ${OPENGL_LIBRARIES})
  target_link_libraries(demo1 ${OPENGL_LIBRARIES})
  target_link_libraries(demo2 ${OPENGL_LIBRARIES})
  target_link_libraries(voronoi ${OPENGL_LIBRARIES})

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
