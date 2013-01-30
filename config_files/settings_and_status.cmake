if(NOT SETTINGS_ADDED EQUAL 1)
  # display status message for important variables
  MESSAGE( STATUS )
  MESSAGE( STATUS "-------------------------------------------------------------------------------" )
  message( STATUS "PROJECT_NAME = ${PROJECT_NAME}")
  MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
  MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
  MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
  message( STATUS "COMPILESPEC_SET = " ${COMPILESPEC_SET})
  if("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    message(STATUS "CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
    #  message(STATUS 
  else()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
      message(STATUS "CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
    endif()  
  endif()
  message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message(STATUS "CMAKE_SHARED_LINKER_FLAGS = ${CMAKE_SHARED_LINKER_FLAGS}")
  MATH(EXPR DISK_BUFFER_SIZE_MB "${DISK_BUFFER_SIZE} / (1024 * 1024)")
  message( STATUS "DISK_BUFFER_SIZE = ${DISK_BUFFER_SIZE_MB} MB")
  message (STATUS "MEMORY_CONSUMPTION_LEVEL = ${MEMORY_CONSUMPTION_LEVEL}")
  MESSAGE( STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}" )
  MESSAGE( STATUS "CMPH_LIBS = ${CMPH_LIBS}")
  message(STATUS "TBB_LIBRARY = ${TBB_LIBRARY}")
  message(STATUS "TBB_INCLUDES = ${TBB_INCLUDES}")
#  message( STATUS "TURBO_LIBS = ${TURBO_LIBS}")
#  message( STATUS "TURBO_INCLUDES = ${TURBO_INCLUDES}")
  message( STATUS "PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
  message( STATUS "CMAKE_SYSTEM_VERSION = ${CMAKE_SYSTEM_VERSION}")
  message( STATUS "CMAKE_SYSTEM_PROCESSOR = ${CMAKE_SYSTEM_PROCESSOR}")
  message( STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
  message( STATUS "CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
  MESSAGE( STATUS "${PROJECT_NAME}_DEPENDS = \"${${PROJECT_NAME}_DEPENDS}\"" )
  MESSAGE( STATUS "BUILD_WITH = \"${BUILD_WITH}\"" )
  MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
  MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
  MESSAGE( STATUS "-------------------------------------------------------------------------------" )
  MESSAGE( STATUS )

  # force some variables that could be defined in the command line to be written to cache
  SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
  SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install ${PROJECT_NAME}" FORCE )
  SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
  SET( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" CACHE PATH
    "Path to custom CMake Modules" FORCE )
  SET( INSTALL_DOC "${INSTALL_DOC}" CACHE BOOL
    "Set to OFF to skip build/install Documentation" FORCE )

  # export build settings
  INCLUDE(CMakeExportBuildSettings)
  CMAKE_EXPORT_BUILD_SETTINGS( "${PROJECT_NAME}BuildSettings.cmake" )
  # export library dependencies (keep this as the last line in the file)
  EXPORT_LIBRARY_DEPENDENCIES( "${PROJECT_NAME}LibDeps.cmake" )
  set(SETTINGS_ADDED 1)
endif()