# To compile the code, make a seperate build folder and then do the following:
#
# cd /location/to/build
# cmake /location/to/symphony/src
# make
#
# If a different compiler than the one cmake detects is needed, then the
# recommended way is to do the following:
#
# cmake -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpic++ /location/to/symphony/src 

cmake_minimum_required(VERSION 2.8)

project(symphony)

# ------------------------------USER OPTIONS----------------------------------#
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -DDEBUG")

# Custon GSL install directory
#set(GSL_ROOT_DIR "")
# ------------------------------END OF USER OPTIONS---------------------------#
# Needed to link python interfaces to precompiled .so libraries
set(BUILD_SHARED_LIBS "ON")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Need to find all packages before starting compilation, so that all headers
# needed by the code are included.

# Should the external libraries needed by the internal libraries go into their
# own CMakeLists in their own folders?

# TODO: Only find packages depending on what options are set.

# System libraries
find_library(MATH_LIBRARIES m REQUIRED)

# External packages
find_package(GSL REQUIRED)
find_package(CBLAS REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)

include(UseCython)

include_directories(${GSL_INCLUDE_DIRS}  ${PYTHON_NUMPY_INCLUDE_DIR} .)

add_subdirectory(maxwell_juettner)
add_subdirectory(power_law)
add_subdirectory(kappa)

add_executable(symphony symphony.c symphony.h bessel_mod.c  calc.c  calc.h
               fits.c  fits.h  integrate.c  integrate.h params.c  params.h)

target_link_libraries(symphony maxwell_juettner power_law kappa
                      ${MATH_LIBRARIES}
                      ${GSL_LIBRARIES}
                      ${CBLAS_LIBRARIES})

cython_add_module(symphonyPy symphonyPy.pyx calc.c integrate.c bessel_mod.c
  params.c params.h fits.c fits.h)
target_link_libraries(symphonyPy kappa power_law maxwell_juettner
  ${GSL_LIBRARIES} ${CBLAS_LIBRARIES})

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("C Compiler       : " ${CMAKE_C_COMPILER})
message("C_FLAGS          : " ${CMAKE_C_FLAGS})
message("GSL Libraries    : " ${GSL_LIBRARIES})
message("CBLAS Libraries  : " ${CBLAS_LIBRARIES})
message("NumPy dir        : " ${PYTHON_NUMPY_INCLUDE_DIR})
message("Cython           : " ${CYTHON_EXECUTABLE})
message("")
