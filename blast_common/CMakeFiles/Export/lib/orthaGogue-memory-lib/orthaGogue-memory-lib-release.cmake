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
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/usr/lib/libtbb.so;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/terminal_input/libcmd_line.a;/usr/local/lib/libcmph.so;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/../log_builder/liblog_builder.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthaGogue/libblast_memory.a"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS blast_memory )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_blast_memory "${_IMPORT_PREFIX}/lib/orthaGogue/libblast_memory.a" )

# Import target "orthaGogue-memory-lib" for configuration "Release"
SET_PROPERTY(TARGET orthaGogue-memory-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(orthaGogue-memory-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "blast_memory;/usr/lib/libtbb.so;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/terminal_input/libcmd_line.a;/usr/local/lib/libcmph.so;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/../log_builder/liblog_builder.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthaGogue/liborthaGogue-memory-lib.a"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS orthaGogue-memory-lib )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_orthaGogue-memory-lib "${_IMPORT_PREFIX}/lib/orthaGogue/liborthaGogue-memory-lib.a" )

# Loop over all imported files and verify that they actually exist
FOREACH(target ${_IMPORT_CHECK_TARGETS} )
  FOREACH(file ${_IMPORT_CHECK_FILES_FOR_${target}} )
    IF(NOT EXISTS "${file}" )
      MESSAGE(FATAL_ERROR "The imported target \"${target}\" references the file
   \"${file}\"
but this file does not exist.  Possible reasons include:
* The file was deleted, renamed, or moved to another location.
* An install or uninstall procedure did not complete successfully.
* The installation package was faulty and contained
   \"${CMAKE_CURRENT_LIST_FILE}\"
but not all the files it references.
")
    ENDIF()
  ENDFOREACH()
  UNSET(_IMPORT_CHECK_FILES_FOR_${target})
ENDFOREACH()
UNSET(_IMPORT_CHECK_TARGETS)

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)