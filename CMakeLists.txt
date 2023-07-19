cmake_minimum_required(VERSION 3.5)

project(psfex VERSION "3.0.0" LANGUAGES C)

set(DEFAULT_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include(GNUInstallDirs)

include(FindPkgConfig)
pkg_search_module(FFTW3 REQUIRED fftw3f)

find_package(GSL REQUIRED)

set(SOURCE_FILES
    src/levmar/Axb.c
    src/levmar/lm.c
    src/levmar/lmbc.c
    src/levmar/lmblec.c
    src/levmar/lmbleic.c
    src/levmar/lmlec.c
    src/levmar/misc.c
    src/wcs/poly.c
    src/context.c
    src/diagnostic.c
    src/dummies.c
    src/fft.c
    src/field_utils.c
    src/field.c
    src/fits/fitsmisc.c
    src/homo.c
    src/lapack_stub.c
    src/makeit2.c
    src/misc.c
    src/pca.c
    src/prefs.c
    src/psf.c
    src/sample_utils.c
    src/vignet.c
    src/wcs_utils.c
    src/xml.c
)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:src>
    PRIVATE
        src
        ${FFTW3_INCLUDE})

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        GSL::gslcblas
        fftw3 fftw3f)

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)
