cmake_minimum_required (VERSION 2.8.2)
project(main)

#MPI
find_package(MPI)

IF(NOT MPI_CXX_FOUND)
	set(MPI_CXX_COMPILER $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/bin/mpicxx)
	set(MPI_C_COMPILER $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/bin/mpicc)
	set(MPI_FORTRAN_COMPILER $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/bin/mpifort)
	find_package(MPI REQUIRED)
ENDIF(NOT MPI_CXX_FOUND)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")

#Trilinos
FIND_PACKAGE(Trilinos PATHS $ENV{TRILINOS_DIR}/lib/cmake/Trilinos $ENV{TRILINOS_DIR})

IF(Trilinos_FOUND)
	MESSAGE("\nFound Trilinos\n")
ENDIF()

#PETSc
find_library(PETSC_LIBRARY 
						NAMES petsc
						HINTS "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib")

#PetIGA
find_library(PETIGA_LIBRARY 
						NAMES petiga
						HINTS "$ENV{PETIGA_DIR}/$ENV{PETSC_ARCH}/lib")

#Sacado
#find_library(SACADO_LIBRARY 
#						NAMES sacado
#						HINTS "$ENV{TRILINOS_DIR}/lib/")


#Include directories
include_directories("${MPI_CXX_INCLUDE_PATH}" 
			"$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include" 
			"$ENV{PETSC_DIR}/include" 
			"$ENV{PETIGA_DIR}/$ENV{PETSC_ARCH}/include" 
			"$ENV{PETIGA_DIR}/include" 
			${Trilinos_TPL_INCLUDE_DIRS}
			${Trilinos_INCLUDE_DIRS}
			"../../../include")

#Add src
file(GLOB SOURCES "../../../src/*.cc" "*.cc")
add_executable(main ${SOURCES})

#Link libraries
LINK_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS} ${Trilinos_INCLUDE_DIRS})
target_link_libraries(main ${PETIGA_LIBRARY} ${PETSC_LIBRARY} ${MPI_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${Trilinos_LIBRARIES})#${SACADO_LIBRARY})

#Debug or release
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(debugRelease
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=RelWithDebInfo ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
