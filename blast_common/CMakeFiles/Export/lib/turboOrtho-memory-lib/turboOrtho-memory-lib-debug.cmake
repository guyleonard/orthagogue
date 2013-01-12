#----------------------------------------------------------------
# Generated CMake target import file for configuration "DEBUG".
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Compute the installation prefix relative to this file.
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)

# Import target "blast_memory" for configuration "DEBUG"
SET_PROPERTY(TARGET blast_memory APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
SET_TARGET_PROPERTIES(blast_memory PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "tbb;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/terminal_input/libcmd_line.a;/usr/local/lib/libcmph.so;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/orthaGogue/libblast_memory.a"
  )

# Import target "orthaGogue-memory-lib" for configuration "DEBUG"
SET_PROPERTY(TARGET orthaGogue-memory-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
SET_TARGET_PROPERTIES(orthaGogue-memory-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "blast_memory;tbb;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/terminal_input/libcmd_line.a;/usr/local/lib/libcmph.so;/home/klatremus/Dokumenter/Work/orthaGogue/orthaGogue_0.9.9.0/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/orthaGogue/liborthaGogue-memory-lib.a"
  )

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
