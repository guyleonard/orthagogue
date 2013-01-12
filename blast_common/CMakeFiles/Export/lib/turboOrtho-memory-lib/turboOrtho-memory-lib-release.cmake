#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Compute the installation prefix relative to this file.
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)

# Import target "blast_memory" for configuration "Release"
SET_PROPERTY(TARGET blast_memory APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(blast_memory PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "tbb;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/terminal_input/libcmd_argument.a;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/terminal_input/libcmd_list.a;/usr/local/lib/libcmph.so;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthaGogue/libblast_memory.a"
  )

# Import target "orthaGogue-memory-lib" for configuration "Release"
SET_PROPERTY(TARGET orthaGogue-memory-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(orthaGogue-memory-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "blast_memory;tbb;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/terminal_input/libcmd_argument.a;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/terminal_input/libcmd_list.a;/usr/local/lib/libcmph.so;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthaGogue/liborthaGogue-memory-lib.a"
  )

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
