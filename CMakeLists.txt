##-----------------------------------------------------------
##
##    Copyright (C) 2015 by the deal.II authors
##
##    This file is part of the deal.II library.
##
##    The deal.II library is free software; you can use it, redistribute
##    it, and/or modify it under the terms of the GNU Lesser General
##    Public License as published by the Free Software Foundation; either
##    version 2.1 of the License, or (at your option) any later version.
##    The full text of the license can be found in the file LICENSE at
##    the top level of the deal.II distribution.
##
##-----------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

##
# CMake script for large project, with one executable
# and a library of objects, depending on deal.II.
#
# The structure of the directory is assumed to be: 
# ./source: containing implementations and one file with "int main()"
# ./include: containing all class declarations
# ./tests: containing googletest
##

# Set the name of the project and target
# If your application follows the structure above, you don't need to 
# specify anything else. 
SET(TARGET ndim_subdivision)
SET(_main source/main.cc)

# Set the _main variable to empty if you don't want an executable
# but only a library
#SET(_main "")

############################################################
# Normally you shouldn't need to change anything below.
############################################################
# Declare all source files the target consists of:
file(GLOB _files source/*cc include/*.h)

# Don't compile the main file into the library
IF(NOT "${_main}" STREQUAL "") 
  LIST(REMOVE_ITEM _files
    ${CMAKE_CURRENT_SOURCE_DIR}/${_main}
    )
ENDIF()

INCLUDE_DIRECTORIES(include)

FIND_PACKAGE(deal.II 9.1.0 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(${TARGET})

OPTION(ENABLE_APP_TESTING "Enable google test in this application." ON)

# We add one library and one target for each type of deal.II library
# we found. If you compiled deal.II with both Release and Debug
# mode, this will generate both Release and Debug programs for you.
# The debug library and program are postfixed with ".g"
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")

set(CMAKE_MACOSX_RPATH 0)

add_definitions(-DRESOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/resources\")

FOREACH(_build_type ${_d2_build_types})
  # Postfix to use everywhere
  SET(_p "${${_build_type}_postfix}")
  # Only build this type, if deal.II was compiled with it.
  IF(CMAKE_BUILD_TYPE MATCHES "${_build_type}" AND
      DEAL_II_BUILD_TYPE MATCHES "${_build_type}")

    MESSAGE("-- Found ${_build_type} version of deal.II.")

    SET(_lib "${TARGET}-lib${_p}")
    MESSAGE("-- Configuring library ${_lib}")

    STRING(TOUPPER "${_build_type}" _BUILD_TYPE)

    ADD_LIBRARY(${_lib} SHARED ${_files})
    FOREACH(_dim 1 2 3)
	SET(_exe "${TARGET}_${_dim}d${${_build_type}_postfix}")
	MESSAGE("-- Configuring executable ${_exe}")
	ADD_EXECUTABLE(${_exe} ${_main})
	TARGET_COMPILE_DEFINITIONS(${_exe} PUBLIC DEAL_II_DIMENSION=${_dim})
	
	TARGET_LINK_LIBRARIES(${_exe} ${_lib})
	DEAL_II_SETUP_TARGET(${_exe} ${_BUILD_TYPE})
    ENDFOREACH()

    SET(TEST_LIBRARIES_${_BUILD_TYPE} ${_lib})

    DEAL_II_SETUP_TARGET(${_lib} ${_BUILD_TYPE})
  ENDIF()
ENDFOREACH()


if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests AND ${ENABLE_APP_TESTING})
    enable_testing()
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIR})
    ADD_SUBDIRECTORY(tests)
endif()


# Doxygen
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
    endif(DOXYGEN_FOUND)
endif()

add_custom_target(indent ./scripts/indent WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# additional target to perform clang-format run, requires clang-format
