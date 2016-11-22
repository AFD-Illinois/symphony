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

cmake_minimum_required(VERSION 3.2)

project(symphony)
set(LIBRARY_VERSION 0.1)

# ------------------------------USER OPTIONS----------------------------------#
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3 -g")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -DDEBUG")

# Custom GSL install directory
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
find_package(PythonInterp)

# Derive the Python site-packages directory.
# http://stackoverflow.com/questions/1242904/
set(_cmd
  "from distutils.sysconfig import get_python_lib"
  "print(get_python_lib(plat_specific=True, prefix=''))")
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" "-c" "${_cmd}"
  OUTPUT_VARIABLE PYTHON_SITE_DIR
  ERROR_VARIABLE _pyerr
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(_pyerr)
  message(FATAL_ERROR "Python command failed:\n${_pyerr}")
endif(_pyerr)

# External packages
find_package(GSL REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)

include(UseCython)

include_directories(${GSL_INCLUDES}  ${PYTHON_NUMPY_INCLUDE_DIR} .)

add_library(symphony
bessel_mod.c
distribution_function_common_routines.c
distribution_function_common_routines.h
integrator/integrands.c
integrator/integrands.h
integrator/integrate.c
integrator/integrate.h
fits.c
fits.h
kappa/kappa.c
kappa/kappa.h
kappa/kappa_fits.c
maxwell_juettner/maxwell_juettner.c
maxwell_juettner/maxwell_juettner.h
maxwell_juettner/maxwell_juettner_fits.c
params.c
params.h
power_law/power_law.c
power_law/power_law.h
power_law/power_law_fits.c
symphony.c
symphony.h
)

target_link_libraries(symphony
                      ${MATH_LIBRARIES}
                      ${GSL_LIBRARIES}
                      ${CBLAS_LIBRARIES})

add_executable(demo demo.c)
target_link_libraries(demo symphony)

cython_add_module(symphonyPy symphonyPy.pyx)
target_link_libraries(symphonyPy symphony
  ${GSL_LIBRARIES} ${CBLAS_LIBRARIES})

set_target_properties(
  symphony
  PROPERTIES
    VERSION ${LIBRARY_VERSION}
    SOVERSION ${LIBRARY_VERSION}
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS symphony
  LIBRARY DESTINATION lib)

set_target_properties(
  symphonyPy
  PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS symphonyPy
  LIBRARY DESTINATION ${PYTHON_SITE_DIR})

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("C Compiler          : " ${CMAKE_C_COMPILER})
message("C_FLAGS             : " ${CMAKE_C_FLAGS})
message("GSL Libraries       : " ${GSL_LIBRARIES})
message("CBLAS Libraries     : " ${CBLAS_LIBRARIES})
message("NumPy dir           : " ${PYTHON_NUMPY_INCLUDE_DIR})
message("Cython              : " ${CYTHON_EXECUTABLE})
message("Library install dir : " ${CMAKE_INSTALL_PREFIX})
message("Python install dir  : " ${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE_DIR})
message("")
